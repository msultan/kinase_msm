#!/bin/env python
import os
from .convert_project import extract_project_wrapper, _sanity_tests
from .data_loader import load_yaml_file

def convert_series(yaml_file, ip_view, protein_list = None):
    """
    :param yaml_file: The yaml file to work with
    :param ip_view: ipython view(required)
    :param protein_list: list of proteins, if None then all
    the proteins in yaml_file["protein_list"] are processed
    :return: converted and concatenated trajectories in
    yaml_file["base_dir"]+protein_name+trajectories
    and the stripped files in
    yaml_file["base_dir"]+protein_name+protein_traj
    """

    yaml_file = load_yaml_file(yaml_file)

    if protein_list is None:
        protein_list = yaml_file["protein_list"]

    for protein in protein_list:
        for proj in yaml_file["project_dict"][protein]:

            protein_folder = os.path.join(yaml_file["base_dir"], protein)
            proj_folder = os.path.join(protein_folder, proj)
            top_folder = os.path.join(proj_folder,"topologies")

            _sanity_tests(protein_folder, proj_folder, top_folder)

            extract_project_wrapper(yaml_file, protein,
                                    proj, ip_view)

    return