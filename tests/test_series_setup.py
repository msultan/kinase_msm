#!/bin/env/python
from kinase_msm.series_setup import setup_series_analysis
import os
import yaml
import glob
import time
from mdtraj.utils.contextmanagers import enter_temp_directory

def create_fake_series():
    os.mkdir("fake_series")
    os.mkdir("./fake_series/fake_kinase1")
    os.mkdir("./fake_series/fake_kinase2")
    os.mkdir("./fake_series/fake_kinase1/fake_proj1")
    os.mkdir("./fake_series/fake_kinase1/fake_proj2")
    os.mkdir("./fake_series/fake_kinase2/fake_proj3")
    return

def test_setup_series_analysis():

    base_dir = os.path.join("./fake_series")
    mdl_dir = os.path.join(base_dir,"new_mdl_dir")
    feature_dir = "feature_dir"
    series_name = "fake_series"
    protein_list = ["fake_kinase1", "fake_kinase2"]
    project_dict = {"fake_kinase1": ["fake_proj1", "fake_proj2"],
                    "fake_kinase2": ["fake_proj3"]}
    mdl_params = {'tica__n_components': 1, 'tica__lag_time': 2,
                  'tica__kinetic_mapping': True, 'tica__shrinkage': 0.01,
                  'cluster__n_clusters': 174}

    with enter_temp_directory():
        create_fake_series()
        setup_series_analysis(base_dir, mdl_dir, feature_dir,
                              series_name, protein_list,
                              project_dict, mdl_params)

        assert os.path.isdir(mdl_dir)
        for protein in protein_list:
            assert os.path.isdir(os.path.join(mdl_dir, protein))

        assert(os.path.isfile(os.path.join(base_dir,"series.yaml")))
        fin = open(os.path.join(mdl_dir,"project.yaml"), 'r')
        yaml_file = yaml.load(fin)

        assert yaml_file["base_dir"] == base_dir
        assert yaml_file["series_name"] == series_name
        assert yaml_file["protein_list"] == protein_list
        assert yaml_file["project_dict"] == project_dict
        assert yaml_file["mdl_params"] == mdl_params

    return

def test_multiple_mdls():
    base_dir = os.path.join("./fake_series")
    mdl_dir = os.path.join(base_dir,"new_mdl_dir")
    feature_dir = "feature_dir"
    series_name = "fake_series"
    protein_list = ["fake_kinase1", "fake_kinase2"]
    project_dict = {"fake_kinase1": ["fake_proj1", "fake_proj2"],
                    "fake_kinase2": ["fake_proj3"]}
    mdl_params = {'tica__n_components': 4, 'tica__lag_time': 223,
                  'tica__kinetic_mapping': True, 'tica__gamma': 0.0121,
                  'cluster__n_clusters': 212}

    with enter_temp_directory():
        create_fake_series()
        for i in range(3):
            setup_series_analysis(base_dir, mdl_dir, feature_dir,
                                  series_name, protein_list,
                                  project_dict, mdl_params)
            time.sleep(1)
        assert len(glob.glob("./fake_series/*/project.yaml")) == 3
    return

