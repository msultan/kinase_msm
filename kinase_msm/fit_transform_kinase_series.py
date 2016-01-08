#!/bin/env python

from __future__ import print_function
from msmbuilder.decomposition import tICA
from msmbuilder.utils import verboseload, verbosedump
import glob
from msmbuilder.msm import BayesianMarkovStateModel, MarkovStateModel
import os
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.dataset import _keynat as keynat
from kinase_msm.data_loader import change_protein_data_dir, change_protein_mdl_dir, load_yaml_file

def fit_protein_tica(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    mdl_params = yaml_file["mdl_params"]

    tica__lag_time = mdl_params["tica__lag_time"]
    tica__gamma = mdl_params["tica__gamma"]
    tica__weighted_transform = mdl_params["tica__weighted_transform"]
    tica__n_components = mdl_params["tica__n_components"]

    protein_tica_mdl = tICA(n_components=tica__n_components, lag_time=tica__lag_time,
                            gamma=tica__gamma,
                            weighted_transform=tica__weighted_transform)

    for protein in yaml_file["protein_list"]:
        print("Fitting to protein %s" % protein)
        change_protein_data_dir(yaml_file["base_dir"], protein)
        featurized_traj = sorted(glob.glob("./%s/*.jl"%
                                 yaml_file["feature_dir"]), key=keynat)
        for f in featurized_traj:
            featurized_path = verboseload(f)
            try:
                protein_tica_mdl.partial_fit(featurized_path)
            except:
                pass
        print("Done partial fitting to protein %s" % protein)
    #dumping the tica_mdl
    tica_mdl_path = os.path.join(mdl_dir, "tica_mdl.pkl")
    verbosedump(protein_tica_mdl, tica_mdl_path)
    return


def transform_protein_tica(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    tica_obj_path = os.path.join(mdl_dir, "tica_mdl.pkl")
    protein_tica_mdl = verboseload(tica_obj_path)
    for protein in yaml_file["protein_list"]:
        change_protein_data_dir(yaml_file["base_dir"], protein)
        print("Transforming protein %s" % protein)
        featurized_traj = sorted(glob.glob("./%s/*.jl"%
                                 yaml_file["feature_dir"]), key=keynat)
        tica_data = {}
        for f in featurized_traj:
            featurized_path = verboseload(f)
            try:
                tica_data[os.path.basename(f)] = \
                    protein_tica_mdl.partial_transform(featurized_path)
            except:
                pass
        change_protein_mdl_dir(yaml_file["base_dir"], protein)
        verbosedump(tica_data, 'tica_data.pkl')
        print("Done transforming protein %s" % protein)
    return


def fit_protein_kmeans(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    mdl_params = yaml_file["mdl_params"]
    cluster__n_clusters = mdl_params["cluster__n_clusters"]

    kmeans_mdl = MiniBatchKMeans(cluster__n_clusters,
                                 batch_size=max(10000, 10 * cluster__n_clusters))
    data = []

    for protein in yaml_file["protein_list"]:
        change_protein_mdl_dir(yaml_file["base_dir"], protein)
        tica_data = verboseload("tica_data.pkl")
        # get all traj
        sorted_list = sorted(tica_data.keys(), key=keynat)
        data.extend([tica_data[i] for i in sorted_list])

    kmeans_mdl.fit(data)
    kmeans_mdl_path = os.path.join(mdl_dir, "kmeans_mdl.pkl")
    try:
        os.remove(kmeans_mdl_path)
    except:
        pass
    verbosedump(kmeans_mdl, kmeans_mdl_path)
    return


def transform_protein_kmeans(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    kmeans_mdl_path = os.path.join(mdl_dir, "kmeans_mdl.pkl")
    kmeans_mdl = verboseload(kmeans_mdl_path)
    for protein in yaml_file["protein_list"]:
        print("Assigning protein %s" % protein)
        change_protein_mdl_dir(yaml_file["base_dir"],protein)
        tica_data = verboseload("tica_data.pkl")
        # do assignments
        assignments = {}
        for i in tica_data.keys():
            assignments[i] = kmeans_mdl.predict([tica_data[i]])[0]
        verbosedump(assignments, 'assignments.pkl')

        print("Done assigning %s" % protein)
    return


def fit_msms(yaml_file):
    mdl_params = yaml_file["mdl_params"]
    msm__lag_time = mdl_params["msm__lag_time"]
    for protein in yaml_file["protein_list"]:
        print(protein)
        change_protein_mdl_dir(yaml_file["base_dir"], protein)
        assignments = verboseload("assignments.pkl")
        msm_mdl = MarkovStateModel(
            lag_time=msm__lag_time, verbose=True).fit(
            [assignments[i] for i in assignments.keys()])
        verbosedump(msm_mdl, "msm_mdl.pkl")
        fixed_assignments = {}
        for i in assignments.keys():
            fixed_assignments[i] = msm_mdl.transform(
                assignments[i], mode='fill')[0]
        verbosedump(fixed_assignments, 'fixed_assignments.pkl')
    return


def fit_bayes_msms(yaml_file):
    mdl_params = yaml_file["mdl_params"]
    msm__lag_time = mdl_params["msm__lag_time"]
    if "bayesmsm__n_samples"  in mdl_params.keys():
        bayesmsm__n_samples = mdl_params["bayesmsm__n_samples"]
    else:
        bayesmsm__n_samples = 800
    if "bayesmsm__n_steps" in mdl_params.keys():
        bayesmsm__n_steps = mdl_params["bayesmsm__n_steps"]
    else:
        bayesmsm__n_steps = 1000000

    for protein in yaml_file["protein_list"]:
        print(protein)
        change_protein_mdl_dir(yaml_file["base_dir"], protein)
        assignments = verboseload("assignments.pkl")
        msm_mdl = BayesianMarkovStateModel(n_samples=bayesmsm__n_samples,
                                           n_steps=bayesmsm__n_steps,
                                           lag_time=msm__lag_time, verbose=True).fit(
            [assignments[i] for i in assignments.keys()])
        verbosedump(msm_mdl, "bayesmsm_mdl.pkl")
    return


def fit_pipeline(base_dir, mdl_dir=None):
    os.chdir(base_dir)
    if mdl_dir is None:
        mdl_dir = os.path.join(base_dir, "mdl_dir")
    yaml_file = load_yaml_file(os.path.join(mdl_dir, "project.yaml"))

    fit_protein_tica(yaml_file)
    transform_protein_tica(yaml_file)
    fit_protein_kmeans(yaml_file)
    transform_protein_kmeans(yaml_file)
    fit_msms(yaml_file)
    fit_bayes_msms(yaml_file)

    return
