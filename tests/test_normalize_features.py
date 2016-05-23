#!/bin/env/python

from kinase_msm.featurize_project import featurize_project_wrapper
from kinase_msm.normalize_features import normalize_project_series
from kinase_msm.data_loader import load_yaml_file
from test_convert_series import _setup_test, _cleanup_test
from nose.tools import with_setup
from msmbuilder.utils import verboseload
import glob
import os
import numpy as np
import mdtraj as mdt

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))

#test to make sure features are normalized.T


@with_setup(_setup_test, _cleanup_test)
def test_normalize_features():

    print(base_dir)
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    normalize_project_series(yaml_file, stride=1)

    all_data=[]

    for kinase in ["kinase_1","kinase_2"]:
        flist = glob.glob(os.path.join(base_dir, "%s"%kinase , "normalized_features/*.jl"))

        for i in flist:
            all_data.extend(verboseload(i))

    assert(np.alltrue(np.isclose(np.mean(all_data, axis=1), 0 , atol=1e-1)))
    assert(np.alltrue(np.isclose(np.std(all_data, axis=1), 1 , atol=0.2)))


