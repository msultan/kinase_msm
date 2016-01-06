import os
from msmbuilder.utils import verboseload, verbosedump
import mdtraj as mdt
import numpy as np

'''
script to load pertinent data for a given protein
and perform some sanity checks
'''


def change_protein_data_dir(base_dir, protein):
    """
    :param base_dir: The base directory for the project
    :param protein: The protein within the project
    :return: Nothing but the cwd should be the protein data dir
    i.e base_dir/protein
    """
    if os.getcwd() != os.path.join(base_dir, protein):
        os.chdir(os.path.join(base_dir, protein))
    return


def change_protein_mdl_dir(base_dir, protein):
    """
    :param base_dir: The base directory for the project
    :param protein: The protein within the project
    :return: Nothing but the cwd should be the protein mdl dir
    base_dir/mdl_dir/protein
    """
    if os.getcwd() != os.path.join(base_dir, "mdl_dir", protein):
        os.chdir(os.path.join(base_dir, "mdl_dir", protein))
    return


def load_traj(base_dir, protein, filename):
    """
    :param base_dir: Project's base dir
    :param protein: Protein of interest
    :param filename: file to load
    :return: the trajectory obj
    """
    change_protein_data_dir(base_dir, protein)
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
    return load_traj(base_dir,protein,filename)[frame_index]


def sanity_test(base_dir, protein, msm_mdl, tica_data, kmeans_mdl, assignments):
    tics_to_use = kmeans_mdl.cluster_centers_.shape[1]
    for i, v in enumerate(tica_data.keys()[:20]):
        # skip
        if not np.isnan(assignments[v]).any():

            assert((msm_mdl.transform(kmeans_mdl.transform(
                [tica_data[v][:, :tics_to_use]])) == assignments[v]).all())
            trj = load_traj(base_dir, mutant, v.split(".h5")[0])
            assert(trj.n_frames == assignments[v].shape[0])
    return


def load_current_protein_model(base_dir, protein, sanity=True):
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
    mdl_dir = os.path.join(base_dir, "mdl_dir")
    prot_mdl_dir = os.path.join(base_dir, protein)

    #load the project level information first
    kmeans_mdl = verboseload(os.path.join(mdl_dir, "kmeans_mdl.pkl"))
    tica_mdl = verboseload(os.path.join(mdl_dir, "tica_mdl.pkl"))

    #now load the protein level information
    tica_data = verboseload(os.path.join(prot_mdl_dir, "tica_data.pkl"))
    # need the fixed assignments because otherwise we will have issues
    assignments = verboseload(os.path.join(prot_mdl_dir, "fixed_assignments.pkl"))
    msm_mdl = verboseload(os.path.join(prot_mdl_dir, "msm_mdl.pkl"))
    # some sanity tests
    if sanity:
        sanity_test(base_dir, protein, msm_mdl,
                    tica_data, kmeans_mdl, assignments)
    return base_dir, mdl_dir, msm_mdl, tica_mdl, tica_data, kmeans_mdl, assignments
