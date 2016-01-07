#!/bin/env python

from __future__ import print_function
import os
import glob
from kinase_msm.series_setup import setup_series_analysis
from kinase_msm.fit_transform_kinase_series import *

def test_pipeline():
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
    mdl_dir = os.path.join(base_dir,"mdl_dir")
    feature_dir = "feature_dir"
    series_name = "fake_series"
    kinase_list = ["kinase_1", "kinase_2"]
    project_dict = {"kinase_1": ["fake_proj1",],
                    "kinase_2": ["fake_proj2"]}
    mdl_params = {'tica__n_components': 1, 'tica__lag_time': 1,
                  'tica__weighted_transform': True, 'tica__gamma': 0.01,
                  'cluster__n_clusters': 2,'msm__lag_time': 1}

    yaml_file = setup_series_analysis(base_dir, mdl_dir, feature_dir,
                              series_name, kinase_list,
                              project_dict, mdl_params)
    fit_protein_tica(yaml_file)
    transform_protein_tica(yaml_file)
    fit_protein_kmeans(yaml_file)
    transform_protein_kmeans(yaml_file)
    fit_msms(yaml_file)
    
    raw_count_obs = 0 
    for p in kinase_list:
        for j in glob.glob(os.path.join(base_dir,p,feature_dir,"*.jl")):
            raw_count_obs += verboseload(j).shape[0]
    tica_mdl = verboseload(os.path.join(mdl_dir,"tica_mdl.pkl"))
    #make sure the mdl is seeing all the data, could probably have a far stronger test here
    assert tica_mdl.n_observations_ == raw_count_obs 
    assert os.path.exists(os.path.join(mdl_dir,"kinase_1/tica_data.pkl"))
    assert os.path.exists(os.path.join(mdl_dir,"kinase_2/tica_data.pkl"))
    assert os.path.exists(os.path.join(mdl_dir,"kinase_1/msm_mdl.pkl"))
    assert os.path.exists(os.path.join(mdl_dir,"kinase_2/msm_mdl.pkl"))
    assert os.path.exists(os.path.join(mdl_dir,"kmeans_mdl.pkl"))
    
    return 
