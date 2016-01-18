#!/bin/env/python

from kinase_msm.featurize_project import featurize_project_wrapper
from kinase_msm.data_loader import load_yaml_file
from multiprocessing.pool import Pool
from test_convert_series import _setup_test, _cleanup_test
from nose.tools import with_setup
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.utils import verboseload
import glob
import os
import numpy as np
import mdtraj as mdt

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))



@with_setup(_setup_test, _cleanup_test)
def test_dihedral_feat():

    print(base_dir)
    pool = Pool(6)
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))

    for prt in ["kinase_1", "kinase_2"]:
        prj = yaml_file["project_dict"][prt][0]
        featurize_project_wrapper(yaml_file, prt, prj,
                                  feat=None, stride=1, view=pool)

        feat = DihedralFeaturizer(types=['phi', 'psi','chi1'])
        flist = glob.glob(os.path.join(base_dir, prt , "protein_traj/*.hdf5"))
        for i in np.random.choice(flist, 3):
            trj = mdt.load(i)
            my_feat = feat.partial_transform(trj)
            expected_fname = os.path.join(base_dir, prt,
                                          yaml_file["feature_dir"],
                                          os.path.splitext(os.path.basename(i))[0]+".jl")
            calc_feat = verboseload(expected_fname)

            assert np.allclose(my_feat, calc_feat)



    return True

