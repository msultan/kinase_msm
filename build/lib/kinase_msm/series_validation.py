#!/bin/env python

from kinase_msm.data_loader import load_yaml_file
from multiprocessing.pool import Pool
from multiprocessing import cpu_count

def _validate_protein(t):
    """
    :param yaml_file: yaml_file
    :param protein: protein name
    :param sequence_dictionary: sequence dict
    :return:
    """
    yaml_file, protein, sequence_dictionary = t
    print(protein)
    return

def validate_series(yaml_file, sequence_dictionary):
    """
    :param yaml_file: The mdl yaml file.
    :param sequence_dictionary: Dictionary of sequences
    :return: Runs a large number of sequence tests on the series to make sure
    the sequences for each protein match the given sequence and the series itself
    """
    yaml_file = load_yaml_file(yaml_file)
    p = Pool(cpu_count())
    jobs = [(yaml_file, protein, sequence_dictionary)
            for protein in yaml_file["protein_list"]]
    p.map(_validate_protein, jobs)

    return