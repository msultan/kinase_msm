#!/bin/env python
from __future__ import print_function, division
import os
import glob
import sys
import tarfile
from msmbuilder.dataset import _keynat as keynat
from mdtraj.formats.hdf5 import HDF5TrajectoryFile
from mdtraj.utils import six
import mdtraj as md
from mdtraj.utils.contextmanagers import enter_temp_directory

def hdf5_concatenate(job_tuple):
    """Concatenate tar bzipped XTC files created by Folding@Home Core17.
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

    proj_folder, top_folder, db_name, run, clone = job_tuple
    path = os.path.join(proj_folder,"RUN%d/CLONE%d/"%(run,clone))
    top = md.load(os.path.join(top_folder,"%d.pdb"%run))
    output_filename =  os.path.join(proj_folder,"trajectories/%d_%d.hdf5"%(run,clone))

    glob_input = os.path.join(path, "results-*.tar.bz2")
    filenames = glob.glob(glob_input)
    filenames = sorted(filenames, key=keynat)

    if len(filenames) <= 0:
        return

    trj_file = HDF5TrajectoryFile(output_filename, mode='a')

    try:
        trj_file._create_earray(where='/', name='processed_filenames',
                                atom=trj_file.tables.StringAtom(1024), shape=(0,))
        trj_file.topology = top.topology
    except trj_file.tables.NodeError:
        pass

    for filename in filenames:
        if six.b(filename) in trj_file._handle.root.processed_filenames:
        # On Py3, the pytables list of filenames has type byte (e.g. b"hey"),
        # so we need to deal with this via six.
            print("Already processed %s" % filename)
            continue
        with enter_temp_directory():
            print("Processing %s" % filename)
            archive = tarfile.open(filename, mode='r:bz2')
            archive.extract("positions.xtc")
            trj = md.load("positions.xtc", top=top)

            for frame in trj:
                trj_file.write(coordinates=frame.xyz,
                               cell_lengths=frame.unitcell_lengths,
                               cell_angles=frame.unitcell_angles)

            trj_file._handle.root.processed_filenames.append([filename])


def extract_project_wrapper(proj_folder, top_folder, view):
    runs=len(glob.glob(proj_folder+"/RUN*"))
    clones=len(glob.glob(proj_folder+"/RUN0/CLONE*"))

    print("Found %d runs and %d clones in %s"%(runs,clones,proj_folder))
    print("Using %d cores to parallelize"%len(view))

    jobs = [(proj_folder,top_folder,db_name,run,clone)
            for run in range(runs)
            for clone in range(clones)]
    result = view.map_sync(hdf5_concatenate,jobs)


    return result
