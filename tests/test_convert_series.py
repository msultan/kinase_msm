#!/bin/env/python

from kinase_msm.series_setup import setup_series_analysis
import os
from nose.tools import with_setup
import shutil

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))

def _setup_test():
    #remove previous mdl_dir
    try:
        shutil.rmtree(os.path.join(base_dir,"mdl_dir"))
    except:
        pass
    setup_series_analysis(base_dir =base_dir,
                          mdl_dir = os.path.abspath(os.path.join(base_dir,"mdl_dir")),
                          feature_dir = "features",
                          series_name="fake_series",
                          protein_list = ["kinase_1", "kinase_2"],
                          project_dict = {"kinase_1":["fake_proj1","fake_proj2"],
                                          "kinase_2":["fake_proj3"]},
                          mdl_params= {"cluster__n_clusters": 2,
                                       "msm__lag_time": 1,
                                       "tica__shrinkage": 0.005,
                                       "tica__lag_time": 1,
                                       "tica__n_components": 2,
                                       "tica__weighted_transform": True}
                         )

def _cleanup_test():
    return

@with_setup(_setup_test,_cleanup_test)
def test_series():
    return True