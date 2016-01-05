#!/bin/env/python
from kinase_msm.series_setup import setup_series_analysis
import os
import yaml


def test_setup_series_analysis():

    base_dir = os.path.join("./fake_series")
    series_name = "fake_series"
    kinase_list = ["fake_kinase1", "fake_kinase2"]
    project_dict = {"fake_kinase1": ["fake_proj1", "fake_proj2"],
                    "fake_kinase2": ["fake_proj3"]}
    mdl_params = {'tica__n_components': 1, 'tica__lag_time': 2,
                  'tica__weighted_transform': True, 'tica__gamma': 0.01,
                  'cluster__n_clusters': 174}

    setup_series_analysis(base_dir, series_name, kinase_list,
                          project_dict, mdl_params)

    assert os.path.isdir(os.path.join(base_dir, "mdl_dir"))
    for kinase in kinase_list:
        assert os.path.isdir(os.path.join(base_dir, "mdl_dir", kinase))

    fin = open(os.path.join(os.path.join(
        base_dir, "mdl_dir"), "project.yaml"), 'r')
    yaml_file = yaml.load(fin)

    assert yaml_file["base_dir"] == base_dir
    assert yaml_file["series_name"] == series_name
    assert yaml_file["kinase_list"] == kinase_list
    assert yaml_file["project_dict"] == project_dict
    assert yaml_file["mdl_params"] == mdl_params

    return
