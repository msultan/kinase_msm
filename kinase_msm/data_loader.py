import os
from msmbuilder.utils import verboseload, verbosedump
import mdtraj as mdt
import numpy as np
import yaml
import glob
from msmbuilder.dataset import _keynat as keynat
import contextlib
import random

'''
script to load pertinent data for a given protein
and perform some sanity checks
'''

@contextlib.contextmanager
def enter_protein_data_dir(yaml_file, protein):
    """Enters the protein's data directory"""
    cwd = os.getcwd()
    os.chdir(os.path.join(yaml_file["base_dir"], protein))
    yield
    os.chdir(cwd)


@contextlib.contextmanager
def enter_protein_mdl_dir(yaml_file, protein):
    """Enters the protein's data directory"""
    cwd = os.getcwd()
    os.chdir(os.path.join(yaml_file["mdl_dir"], protein))
    yield
    os.chdir(cwd)


def load_random_traj(yaml_file, protein):
    yaml_file = load_yaml_file(yaml_file)
    with enter_protein_data_dir(yaml_file, protein):
        traj_folder = os.path.join(os.getcwd(),"protein_traj")
        traj_files = sorted(glob.glob(os.path.join(traj_folder,"*.hdf5" )),
                        key=keynat)
        trj = load_traj(yaml_file["base_dir"], protein,
                        os.path.basename(random.choice(traj_files)))
    #need just one
    return trj

def load_traj(base_dir, protein, filename):
    """
    :param base_dir: Project's base dir
    :param protein: Protein of interest
    :param filename: file to load
    :return: the trajectory obj
    """
    os.chdir(os.path.join(base_dir, protein))
    filename = os.path.splitext(filename)[0]
    return mdt.load("./protein_traj/%s.hdf5" % filename)


def load_frame(base_dir, protein, filename, frame_index):
    """
    :param base_dir: Project's base dir
    :param protein: Protein of interest
    :param filename: file to load
    :param frame_index: needed frame
    :return: The required frame
    """
    return load_traj(base_dir, protein, filename)[frame_index]


def _sanity_test(base_dir, protein, msm_mdl, tica_data, kmeans_mdl, assignments):
    tics_to_use = kmeans_mdl.cluster_centers_.shape[1]
    for i, v in enumerate(list(tica_data.keys())[:20]):
        # skip
        if not np.isnan(assignments[v]).any():
            assert((msm_mdl.transform(kmeans_mdl.transform(
                [tica_data[v][:, :tics_to_use]])) == assignments[v]).all())
    return


def load_yaml_file(yaml_file):
    if isinstance(yaml_file, dict):
        return yaml_file
    else:
        return yaml.load(open(yaml_file, 'r'))


def load_current_protein_model(yaml_file, protein, sanity=True):
    """
    :param base_dir: Base directory for the project
    :param protein: Protein for which to load
    :param sanity: Whether or not to run sanity tests
    :return: base_dir, mdl_dir,
                msm_mdl, tica_mdl,
                tica_data, kmeans_mdl,
                fixed_assignments for the model currently stored in
                mdl_dir and mdl_dir/protein
    """
    yaml_file = load_yaml_file(yaml_file)
    base_dir = yaml_file["base_dir"]
    mdl_dir = yaml_file["mdl_dir"]

    prot_mdl_dir = os.path.join(mdl_dir, protein)

    # load the project level information first
    kmeans_mdl = verboseload(os.path.join(mdl_dir, "kmeans_mdl.pkl"))
    tica_mdl = verboseload(os.path.join(mdl_dir, "tica_mdl.pkl"))

    # now load the protein level information
    tica_data = verboseload(os.path.join(prot_mdl_dir, "tica_data.pkl"))
    # need the fixed assignments because otherwise we will have issues
    assignments = verboseload(os.path.join(
        prot_mdl_dir, "fixed_assignments.pkl"))
    msm_mdl = verboseload(os.path.join(prot_mdl_dir, "msm_mdl.pkl"))
    # some sanity tests
    if sanity:
        _sanity_test(base_dir, protein, msm_mdl,
                     tica_data, kmeans_mdl, assignments)
    return base_dir, mdl_dir, msm_mdl, tica_mdl, tica_data, kmeans_mdl, assignments
