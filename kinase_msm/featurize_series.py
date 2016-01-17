#!/bin/env python
from kinase_msm.featurize_project import featurize_project_wrapper
from kinase_msm.data_loader import load_yaml_file

def featurize_series(yaml_file, ip_view, protein_list = None):
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
            featurize_project_wrapper(yaml_file, protein,proj, None, 1, ip_view)
    return