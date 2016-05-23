#!/bin/env python
import os
import glob
import mdtraj as mdt
from kinase_msm.data_loader import load_yaml_file
from msmbuilder.dataset import _keynat as keynat
from msmbuilder.utils import verbosedump, verboseload
import pandas as pd
from sklearn import preprocessing
from .featurize_project import _check_output_folder_exists
from .data_loader import enter_protein_mdl_dir, enter_protein_data_dir

#normalize


def normalize_project_series(yaml_file, output_folder="normalized_features",
                             stride=40,nrm=None):
    """
    routine to take a set of proteins features stored in the feature_dir and
    normalize them by removing the mean and setting variance to 1 using the standard
    scaler. The normalizer is dumped into the mdl dir.
    :param yaml_file: The yaml file to work with.
    :param output_folder: The name of the output folder to dump normalized features in
    :param stride: The initial stride in files to fit the normalizer with.
    This is necessary to prevent memory errors. defaults to every 40th file
    :param nrm: previously fit normalizer. else it uses the standard scaler from
    scikitlearn
    :return:
    """
    yaml_file = load_yaml_file(yaml_file)
    #setup normalizer
    if nrm is None:
        nrm = preprocessing.StandardScaler()

    all_data = {}
    for prt in yaml_file["protein_list"]:
        with enter_protein_data_dir(yaml_file, prt):
            print(prt)
            flist = glob.glob("./%s/*.jl"%(yaml_file["feature_dir"]))[::stride]
            for f in flist:
                 all_data[f]=verboseload(f)

    seq=[]
    for i in all_data.keys():
       seq.extend(all_data[i])

    #fit it
    nrm.fit(seq)
    #dump it into the mdl dir.
    verbosedump(nrm,"%s/nrm.h5"%yaml_file["mdl_dir"])

    for prt in yaml_file["protein_list"]:
        _check_output_folder_exists(yaml_file, prt, output_folder)

        with enter_protein_data_dir(yaml_file, prt):
            output_folder = os.path.abspath(output_folder)
            flist = glob.glob("./%s/*.jl"%(yaml_file["feature_dir"]))
            for f in flist:
                res = verboseload(f)
                res = nrm.transform(res)
                verbosedump(res,"%s/%s"%(output_folder,os.path.basename(i)))

    return
