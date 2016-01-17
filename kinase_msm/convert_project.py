#!/bin/env python
from __future__ import print_function, division
import os
import glob
import tarfile
from msmbuilder.dataset import _keynat as keynat
from mdtraj.formats.hdf5 import HDF5TrajectoryFile
from mdtraj.utils import six
import mdtraj as md
from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import load_yaml_file

class HDF5TrajectoryFileWrapper():
    def __init__(self,file):
        assert isinstance(file, HDF5TrajectoryFile)
        self.file = file

    def setup(self, prt_top):
        """
        :param hdf5_file: The hdf5 file to use
        :param prt_top: The protein topology
        :return:
        """
        try:
            self.file._create_earray(where='/', name='processed_filenames',
                                    atom=self.file.tables.StringAtom(1024),
                                    shape=(0,))
            self.file.topology = prt_top
        except self.file.tables.NodeError:
            pass
        return

    def validate_filename(self, index, filename, filenames):
        return True

    def check_filename(self,filename):
        """
        checks if a filename exists in a given hdf5 file
        :param filename: the filename
        :param hdf5_file: the hdf5 file
        :return:
        """
        return six.b(filename) in self.file._handle.root.processed_filenames

    def write_file(self,filename,trj):
        for frame in trj:
            self.file.write(coordinates=frame.xyz,
                       cell_lengths=frame.unitcell_lengths,
                       cell_angles=frame.unitcell_angles)
        self.file._handle.\
                 root.processed_filenames.append([filename])
        return


def _sanity_tests(protein_folder, proj_folder, top_folder):
    """
    :param proj_folder: The project folder for a protein
    :param top_folder: The topology folder for a protein
    :return:
    """
    if not os.path.isdir(top_folder):
        sys.exit("Toplogies Folder Doesnt exist.Exiting!")

    if not os.path.isdir(os.path.join(protein_folder,"trajectories")):
        print("Trajectories folder doesnt exist.Creating")
        os.makedirs(os.path.join(protein_folder,"trajectories"))

    if not os.path.isdir(os.path.join(protein_folder,"protein_traj")):
        print("Trajectories folder doesnt exist.Creating")
        os.makedirs(os.path.join(protein_folder,"protein_traj"))

    return


def _traj_loader(filename, top):
    if os.path.isdir(filename):
        return md.load("%s/positions.xtc"%filename, top=top)
    elif filename.endswith(".bz2"):
        archive = tarfile.open(filename, mode='r:bz2')
        archive.extract("positions.xtc")
        return md.load("positions.xtc", top=top)
    else:
        raise Exception("%s is neither a folder nor a tar.bz2 file")
    return


def hdf5_concatenate(job_tuple):
    """Concatenate tar bzipped or nonbized XTC files created by Folding@Home .
    Parameters
    ----------
    path : str
        Path to directory containing "results-*.tar.bz2".  E.g. a single CLONE directory.
    top : mdtraj.Topology
        Topology for system
    output_filename : str
        Filename of output HDF5 file to generate.
    Notes
    -----
    We use HDF5 because it provides an easy way to store the metadata associated
    with which files have already been processed.
    """

    proj, protein_folder, proj_folder, top_folder, run, clone, protein_only = job_tuple

    path = os.path.join(proj_folder,"RUN%d/CLONE%d/"%(run,clone))
    top = md.load(os.path.join(top_folder,"%d.pdb"%run))
    str_top = top.remove_solvent()


    #output path for stripped trajectory
    strip_prot_out_filename = os.path.join(protein_folder,
                                           "protein_traj/%s_%d_%d.hdf5"%(proj,run,clone))
    str_trj_file = HDF5TrajectoryFile(strip_prot_out_filename, mode='a')
    str_trj_file_wrapper = HDF5TrajectoryFileWrapper(str_trj_file)
    str_trj_file_wrapper.setup(str_top.topology)

    if not protein_only:
        #output path for full trajectory
        output_filename =  os.path.join(protein_folder,
                                    "trajectories/%s_%d_%d.hdf5"%(proj,run,clone))
        trj_file = HDF5TrajectoryFile(output_filename, mode='a')
        trj_file_wrapper = HDF5TrajectoryFileWrapper(trj_file)
        trj_file_wrapper.setup(top.topology)

    glob_input = os.path.join(path, "results*")
    filenames = sorted(glob.glob(glob_input), key=keynat)

    if len(filenames) <= 0:
        return

    for index, filename in enumerate(filenames):
        #if we find it in both then no problem we can continue to the next filename
        if (not protein_only and trj_file_wrapper.check_filename(filename)) and \
                str_trj_file_wrapper.check_filename(filename):
            print("Already processed %s" % filename)
            continue
        with enter_temp_directory():
            print("Processing %s" % filename)
            trj = _traj_loader(filename,top)
            if not protein_only and not trj_file_wrapper.check_filename(filename) :
                if trj_file_wrapper.validate_filename(index, filename, filenames):
                    trj_file_wrapper.write_file(filename, trj)
            #now the stripped file
            if not str_trj_file_wrapper.check_filename(filename):
                if str_trj_file_wrapper.validate_filename(index, filename, filenames):
                    trj = trj.remove_solvent()
                    str_trj_file_wrapper.write_file(filename, trj)

    return


def extract_project_wrapper(yaml_file, protein, proj,
                            view, protein_only = False,):

    yaml_file = load_yaml_file(yaml_file)
    base_dir = yaml_file["base_dir"]

    #get the paths
    protein_folder = os.path.join(base_dir,protein)
    proj_folder = os.path.join(protein_folder, proj)
    top_folder = os.path.join(proj_folder, "topologies")

    _sanity_tests(protein_folder, proj_folder, top_folder)

    #get the runs/clones

    runs = len(glob.glob(proj_folder+"/RUN*"))
    clones = len(glob.glob(proj_folder+"/RUN0/CLONE*"))

    print("Found %d runs and %d clones in %s"%(runs,clones,proj_folder))

    jobs = [(proj, protein_folder, proj_folder, top_folder, run, clone, protein_only)
            for run in range(runs)
            for clone in range(clones)]
    result = view.map(hdf5_concatenate,jobs)


    return result
