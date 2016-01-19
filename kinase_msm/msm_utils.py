#!/bin/evn python
from kinase_msm.data_loader import load_yaml_file
import numpy as np
import os
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from kinase_msm.data_loader import load_frame
from kinase_msm.data_transformer import create_assignment_matrix, create_tics_array


def _sample_state():
    """
    Helper function to do the actual sampling
    :return:
    """
    return

def sample_discarded_states(yaml_file, prt_list=None):
    """
    :param yaml_file: The model yaml file to work with
    :param prt_list:
    :return:
    """
    yaml_file = load_yaml_file(yaml_file)
    if prt_list is None:
        prt_list = yaml_file["protein_list"]

    ser = ProteinSeries(yaml_file)
    for prt in prt_list:
        prt_mdl = Protein(ser, prt)


def sample_states(yaml_file, prt):

    return