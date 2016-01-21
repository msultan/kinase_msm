#!/bin/env/python

from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import enter_protein_data_dir
from test_series_setup import create_fake_series
import os

def test_change_protein_data_dir():
    with enter_temp_directory():
        create_fake_series()
        yaml_file ={}
        yaml_file["base_dir"] = "./fake_series"
        protein = "fake_kinase1"
        with enter_protein_data_dir(yaml_file, protein):
            current_folder_path, current_folder_name = os.path.split(os.getcwd())
            assert current_folder_name == "fake_kinase1"
    return
