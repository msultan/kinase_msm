#!/bin/env/python

from kinase_msm.series_setup import setup_series_analysis
from kinase_msm.convert_project import extract_project_wrapper
from kinase_msm.data_loader import load_yaml_file
from multiprocessing.pool import Pool
from test_convert_series import _setup_test, _cleanup_test
from nose.tools import with_setup
from mdtraj.utils.contextmanagers import enter_temp_directory
import mdtraj as mdt
from msmbuilder.dataset import _keynat as keynat
import glob
import os

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))

def _trj_load(file, top):
    if os.path.isdir(file):
        return mdt.load(os.path.join(file,"positions.xtc"),top=top)
    else:
        os.system("tar -xvjf %s"%file)
        return mdt.load("positions.xtc", top=top)
    return

def _load_project_clone(protein, project, run, clone):
    main_dir = base_dir
    with enter_temp_directory():
        top = mdt.load(os.path.join(main_dir, protein,
                                    project, "topologies",
                                    "%d.pdb"%run))
        t = [_trj_load(f,top) for f in
             sorted(glob.glob(os.path.join(main_dir,
                                           protein,
                                           project,
                                           "RUN%d"%run,"CLONE%d"%clone,
                                           "results*")),
                        key=keynat)]
        print("Length of t is :", len(t))
        print(t[0])
        trj = t[0] + t[1:]

    return trj, trj.remove_solvent()

@with_setup(_setup_test, _cleanup_test)
def test_convert_project():

    print(base_dir)
    pool = Pool(2)
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))

    def test_hdf5(protein, p, r, clone):
        trj, stripped_trj = _load_project_clone(protein, p, r, clone)
        trj2 = mdt.load(os.path.join(base_dir, protein, "trajectories/%s_%d_0.hdf5"%(p,r)))
        trj3 = mdt.load(os.path.join(base_dir, protein,"protein_traj/%s_%d_0.hdf5"%(p,r)))

        for i in ["top", "n_atoms", "n_chains", "n_frames", "n_residues"]:
            assert getattr(trj, i) == getattr(trj2,i)

        for i in ["top", "n_atoms", "n_chains", "n_frames", "n_residues"]:
            assert getattr(stripped_trj, i) == getattr(trj3,i)

        assert (trj.xyz==trj2.xyz).all()
        assert (stripped_trj.xyz==trj3.xyz).all()

        return True



    for i in range(3):
        extract_project_wrapper(yaml_file, "kinase_1", "fake_proj1", pool)
        extract_project_wrapper(yaml_file, "kinase_1", "fake_proj2", pool)

        assert test_hdf5("kinase_1", "fake_proj1", 1, 0)
        assert test_hdf5("kinase_1", "fake_proj1", 1, 0)
        assert test_hdf5("kinase_1", "fake_proj2", 0, 0)

        extract_project_wrapper(yaml_file, "kinase_2", "fake_proj3", pool)
        assert test_hdf5("kinase_2", "fake_proj3", 0, 0)


    return True


