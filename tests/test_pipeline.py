#!/bin/env python

from __future__ import print_function
import os
import glob
import numpy as np
from msmbuilder.utils import verbosedump
from kinase_msm.series_setup import setup_series_analysis
from kinase_msm.fit_transform_kinase_series import *
from mdtraj.utils.contextmanagers import enter_temp_directory


def create_fake_data(base_dir, protein_list, project_dict):
    np.random.seed(42)
    for protein in protein_list:
        os.mkdir(protein)
        os.mkdir(os.path.join(protein, "feature_dir"))
        for project in project_dict[protein]:
            os.mkdir(os.path.join(protein,project))
        for i in range(5):
            X = np.random.randn(20, 3)
            verbosedump(X, os.path.join(protein, "feature_dir" ,"%d.jl"%i))
    return

def test_pipeline():
    with enter_temp_directory():
        base_dir = os.path.abspath(os.path.curdir)
        mdl_dir = os.path.join(base_dir,"mdl_dir")
        feature_dir = "feature_dir"
        series_name = "fake_series"
        protein_list = ["kinase_1", "kinase_2"]
        project_dict = {"kinase_1": ["fake_proj1",],
                        "kinase_2": ["fake_proj2"]}
        mdl_params = {'tica__n_components': 1, 'tica__lag_time': 1,
                  'tica__weighted_transform': True, 'tica__gamma': 0.01,
                  'cluster__n_clusters': 2,'msm__lag_time': 1, 'bayesmsm__n_samples':1,
                  'bayesmsm__n_steps':1}

        create_fake_data(base_dir, protein_list, project_dict)


        yaml_file = setup_series_analysis(base_dir, mdl_dir, feature_dir,
                                  series_name, protein_list,
                                  project_dict, mdl_params)
        fit_protein_tica(yaml_file)
        transform_protein_tica(yaml_file)
        fit_protein_kmeans(yaml_file)
        transform_protein_kmeans(yaml_file)
        fit_msms(yaml_file)
        fit_bayes_msms(yaml_file)

        raw_count_obs = 0
        for p in protein_list:
            for j in glob.glob(os.path.join(base_dir,p,feature_dir,"*.jl")):
                raw_count_obs += verboseload(j).shape[0]
        tica_mdl = verboseload(os.path.join(mdl_dir,"tica_mdl.pkl"))
        #make sure the mdl is seeing all the data, could probably have a far stronger test here
        assert tica_mdl.n_observations_ == raw_count_obs
        assert os.path.exists(os.path.join(mdl_dir,"kinase_1/tica_data.pkl"))
        assert os.path.exists(os.path.join(mdl_dir,"kinase_2/tica_data.pkl"))
        assert os.path.exists(os.path.join(mdl_dir,"kinase_1/msm_mdl.pkl"))
        assert os.path.exists(os.path.join(mdl_dir,"kinase_2/msm_mdl.pkl"))
        assert os.path.exists(os.path.join(mdl_dir,"kinase_2/bayesmsm_mdl.pkl"))
        assert os.path.exists(os.path.join(mdl_dir,"kmeans_mdl.pkl"))

        return
