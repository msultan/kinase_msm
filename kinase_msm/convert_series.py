#!/bin/env python
import os
from kinase_msm.convert_project import extract_project_wrapper
from kinase_msm.data_loader import load_yaml_file

def _sanity_tests(proj_folder, top_folder):
    """
    :param proj_folder: The project folder for a protein
    :param top_folder: The topology folder for a protein
    :return:
    """
    yaml_file = load_yaml_file(yaml_file)
    proj_folder = os.path.join(yaml_file["base_dir"],protein)
    top_folder = os.path.join(proj_folder, "topologies")
    if not os.path.isdir(top_folder):
        sys.exit("Toplogies Folder Doesnt exist.Exiting!")

    if not os.path.isdir(os.path.join(proj_folder,"trajectories")):
        print("Trajectories folder doesnt exist.Creating")
        os.makedir(os.path.join(proj_folder,"trajectories"))

    return

def convert_series(yaml_file, ip_view, protein_list = None):
    """
    :param yaml_file: The yaml file to work with
    :param ip_view: ipython view(required)
    :param protein_list: list of proteins, if None then all
    the proteins in yaml_file["protein_list"] are processed
    :return: converted and concatenated trajectories in
    yaml_file["base_dir"]+protein_name+trajectories
    """

    yaml_file = load_yaml_file(yaml_file)

    if protein_list is None:
        protein_list = yaml_file["protein_list"]

    for protein in protein_list:
        for proj in yaml_file["project_dict"][protein]:
            proj_folder = os.path.join(yaml_file["base_dir"], protein,proj)
            top_folder = os.path.join(proj_folder,"topologies")
            _sanity_tests(yaml_file,protein)
            extract_project_wrapper(proj_folder,
                                top_folder,
                                view=ip_view)
