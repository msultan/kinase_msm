#!/bin/env python

from __future__ import print_function
from msmbuilder.decomposition import tICA
from msmbuilder.utils import verboseload, verbosedump
import glob
from msmbuilder.msm import BayesianMarkovStateModel, MarkovStateModel
from msmbuilder.msm.validation import BootStrapMarkovStateModel
import os
from msmbuilder.cluster import MiniBatchKMeans as KMeans
from msmbuilder.dataset import _keynat as keynat
from .data_loader import enter_protein_data_dir, enter_protein_mdl_dir, load_yaml_file


def fit_protein_tica(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    mdl_params = yaml_file["mdl_params"]

    current_mdl_params={}
    for i in mdl_params.keys():
        if i.startswith("tica__"):
            current_mdl_params[i.strip("tica__")] = mdl_params[i]

    protein_tica_mdl = tICA(**current_mdl_params)

    for protein in yaml_file["protein_list"]:
        print("Fitting to protein %s" % protein)
        with enter_protein_data_dir(yaml_file, protein):
            featurized_traj = sorted(glob.glob("./%s/*.jl" %
                                               yaml_file["feature_dir"]), key=keynat)
            for f in featurized_traj:
                featurized_path = verboseload(f)
                try:
                    protein_tica_mdl.partial_fit(featurized_path)
                except:
                    pass
            print("Done partial fitting to protein %s" % protein)
    # dumping the tica_mdl
    tica_mdl_path = os.path.join(mdl_dir, "tica_mdl.pkl")
    verbosedump(protein_tica_mdl, tica_mdl_path)
    return


def transform_protein_tica(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    tica_obj_path = os.path.join(mdl_dir, "tica_mdl.pkl")
    protein_tica_mdl = verboseload(tica_obj_path)
    for protein in yaml_file["protein_list"]:
        with enter_protein_data_dir(yaml_file, protein):
            print("Transforming protein %s" % protein)
            featurized_traj = sorted(glob.glob("./%s/*.jl" %
                                               yaml_file["feature_dir"]), key=keynat)
            tica_data = {}
            for f in featurized_traj:
                featurized_path = verboseload(f)
                try:
                    tica_data[os.path.basename(f)] = \
                        protein_tica_mdl.partial_transform(featurized_path)
                except:
                    pass
            with enter_protein_mdl_dir(yaml_file, protein):
                verbosedump(tica_data, 'tica_data.pkl')
                print("Done transforming protein %s" % protein)
    return


def fit_protein_kmeans(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    mdl_params = yaml_file["mdl_params"]
    cluster__n_clusters = mdl_params["cluster__n_clusters"]

    kmeans_mdl = KMeans(cluster__n_clusters,batch_size = 100*cluster__n_clusters)
    data = []

    for protein in yaml_file["protein_list"]:
        with enter_protein_mdl_dir(yaml_file, protein):
            tica_data = verboseload("tica_data.pkl")
            # get all traj
            sorted_list = sorted(tica_data.keys(), key=keynat)
            data.extend([tica_data[i] for i in sorted_list])

    kmeans_mdl.fit(data)
    kmeans_mdl_path = os.path.join(mdl_dir, "kmeans_mdl.pkl")
    verbosedump(kmeans_mdl, kmeans_mdl_path)
    return


def transform_protein_kmeans(yaml_file):
    mdl_dir = yaml_file["mdl_dir"]
    kmeans_mdl_path = os.path.join(mdl_dir, "kmeans_mdl.pkl")
    kmeans_mdl = verboseload(kmeans_mdl_path)
    for protein in yaml_file["protein_list"]:
        print("Assigning protein %s" % protein)
        with enter_protein_mdl_dir(yaml_file, protein):
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
        with enter_protein_mdl_dir(yaml_file, protein):
            print(protein)
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

def fit_bootstrap(yaml_file,pool=None):
    mdl_params = yaml_file["mdl_params"]
    msm__lag_time = mdl_params["msm__lag_time"]
    if "bootstrap__n_samples" in mdl_params.keys():
        bootstrap__n_samples = mdl_params["bootstrap__n_samples"]
    else:
        bootstrap__n_samples = 100
    for protein in yaml_file["protein_list"]:
        with enter_protein_mdl_dir(yaml_file, protein):
            print(protein)
            assignments = verboseload("assignments.pkl")
            msm_mdl =BootStrapMarkovStateModel(n_samples=bootstrap__n_samples, n_procs=2,
                                               msm_args ={'lag_time': msm__lag_time}
                                               )
            msm_mdl.fit([assignments[i] for i in assignments.keys()], pool=pool)
            verbosedump(msm_mdl, "bootstrap_msm_mdl.pkl")
            verbosedump(msm_mdl.mle_, "msm_mdl.pkl")
            fixed_assignments = {}
            for i in assignments.keys():
                fixed_assignments[i] = msm_mdl.mle_.transform(
                    assignments[i], mode='fill')[0]
            verbosedump(fixed_assignments, 'fixed_assignments.pkl')
    return            

def fit_bayes_msms(yaml_file):
    mdl_params = yaml_file["mdl_params"]
    msm__lag_time = mdl_params["msm__lag_time"]
    if "bayesmsm__n_samples" in mdl_params.keys():
        bayesmsm__n_samples = mdl_params["bayesmsm__n_samples"]
    else:
        bayesmsm__n_samples = 800
    if "bayesmsm__n_steps" in mdl_params.keys():
        bayesmsm__n_steps = mdl_params["bayesmsm__n_steps"]
    else:
        bayesmsm__n_steps = 1000000

    for protein in yaml_file["protein_list"]:
        with enter_protein_mdl_dir(yaml_file, protein):
            print(protein)
            assignments = verboseload("assignments.pkl")
            msm_mdl = BayesianMarkovStateModel(n_samples=bayesmsm__n_samples,
                                               n_steps=bayesmsm__n_steps,
                                               lag_time=msm__lag_time,
                                               ergodic_cutoff=1.0/msm__lag_time,
                                               verbose=True).fit(
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
    fit_bootstrap(yaml_file)

    return
