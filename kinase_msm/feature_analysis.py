'''
Set of scripts to make feature analysis easier.
'''
import glob
import os
import numpy as np
from msmbuilder.utils import load
from kinase_msm.data_loader import load_yaml_file
from kinase_msm.data_loader import load_random_traj, \
    enter_protein_data_dir, enter_protein_mdl_dir


def pull_features(yaml_file, prt, skip=1, feature_indices=None):
    """
    Simple utility to pull certain features from the feature_folder object
    :param prt: Protein model to use
    :param skip: skip for each file(defaults to 1)
    :param feature_indices: which indices to pull
    :return: dictionary keyed on file name with feature values as arrays
    """
    yaml_file = load_yaml_file(yaml_file)
    all_f ={}
    with enter_protein_data_dir(yaml_file, prt.name):
        feature_file_list = glob.glob("./%s/*.jl"%yaml_file["feature_dir"])
        for i in feature_file_list:
            all_f[os.path.basename(i)]=load(i)[:, feature_indices]

    return all_f