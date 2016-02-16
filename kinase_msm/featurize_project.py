#!/bin/env python

import os
import glob
import mdtraj as mdt
from kinase_msm.data_loader import load_yaml_file
from msmbuilder.dataset import _keynat as keynat
from msmbuilder.utils import verbosedump
import pandas as pd
from msmbuilder.featurizer import DihedralFeaturizer

def featurize_file(job_tuple):

    yaml_file, protein, proj, feat, traj_file,stride = job_tuple
    yaml_file = load_yaml_file(yaml_file)

    if feat is None:
        feat = DihedralFeaturizer(types=['phi', 'psi','chi1'])

    _check_output_folder_exists(yaml_file, protein)

    output_folder = os.path.join(yaml_file["base_dir"],
                                 protein,
                                 yaml_file["feature_dir"])

    traj_name = os.path.splitext(os.path.basename(traj_file))[0]
    output_fname = os.path.join(output_folder, traj_name+".jl")

    feat_descriptor = os.path.join(output_folder, "feature_descriptor.h5")

    trj = mdt.load(traj_file)
    features = feat.partial_transform(trj)
    verbosedump(features, output_fname)

    if not os.path.isfile(feat_descriptor) and hasattr(feat, "describe_features"):
        dih_df = pd.DataFrame(feat.describe_features(trj[0]))
        verbosedump(dih_df, feat_descriptor)

    return


def _check_output_folder_exists(yaml_file, protein, folder_name=None):
    yaml_file = load_yaml_file(yaml_file)
    if folder_name is None:
        folder_name= yaml_file["feature_dir"]
    output_folder = os.path.join(yaml_file["base_dir"],
                                  protein,folder_name)

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    return

def featurize_project_wrapper(yaml_file, protein, proj, feat=None, stride=1, view=None):
    """
    Wrapper function for featurizing project.
    :param yaml_file: The yaml file to work with
    :param protein: Protein Name
    :param proj: Project name
    :param feat: Featurization obj. If none, it defaults to
    phi, psi and chi1. Should support a describe_features attribute
    :param view: ipython view or pool view to parallelize over.
    :return:
    """

    yaml_file = load_yaml_file(yaml_file)
    base_dir = yaml_file["base_dir"]

    _check_output_folder_exists(yaml_file, protein)
    #get the paths
    traj_folder = os.path.join(base_dir, protein, "protein_traj")
    traj_files = sorted(glob.glob(os.path.join(traj_folder,"*.hdf5" )),
                        key=keynat)
    print("Found %d files for featurization in %s"
          %(len(traj_files), traj_folder))

    jobs = [(yaml_file, protein, proj, feat, traj_file, stride) for traj_file in traj_files]

    result = view.map(featurize_file, jobs)


    return result