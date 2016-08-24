#!/bin/evn python
from kinase_msm.data_loader import load_yaml_file, enter_protein_mdl_dir
import numpy as np
import os
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from kinase_msm.data_loader import load_frame
from kinase_msm.data_transformer import create_assignment_matrix, create_tics_array
from msmbuilder.utils.nearest import KDTree
"""
Set of helper scripts for sampling tics
"""

def sample_dimension(data, dimension, n_frames, scheme="linear"):
    """Function to sample a dimension of the data
    using one of 3 schemes. All other dimensions are ignored.

    Parameters
    ----------
    data : list of lists
        List of low dimensional data(output of tica)
    dimension : int
        dimension to sample on
    n_frames: int
        Number of frames required
    scheme: string
        One of either "linear", "random" or "edges". Linear
        samples the tic linearly, random samples randomly
        thereby taking approximate free energies into account,
        and edges samples the edges of the tic only.

    Returns
    -------
       list of tuples where first number is the trajectory index and
       second is the frame index
    """
    d_data = [i[:,dimension][:,np.newaxis] for i in data]

    #sort it because all three sampling schemes use it

    all_vals = []
    for i in d_data:
        all_vals.extend(i.flatten())
    all_vals = np.sort(all_vals)

    #get lineraly placed points
    if scheme=="linear":
        max_val = all_vals[-1]
        min_val = all_vals[0]
        spaced_points = np.linspace(min_val, max_val, n_frames)

    elif scheme=="random":
        spaced_points = np.sort(np.random.choice(all_vals, n_frames))

    elif scheme=="edge":
        _cut_point = np.int(n_frames / 2)
        spaced_points = np.hstack((all_vals[:_cut_point], all_vals[-_cut_point:]))
    else:
        raise ValueError("Scheme has be to one of linear, random or edge")

    tree = KDTree(d_data)

    return_vec = []
    for pt in spaced_points:
        dis, ind = tree.query([pt])
        return_vec.append(ind)

    return return_vec

def sample_region(data, pt_dict, n_frames,):
    """Function to sample a region of the data.

    Parameters
    ----------
    data : list of lists
        List of low dimensional data(output of tica)
    pt_dict : dict
        Dictionary where the keys are the dimensions and the
        value is the value of the dimension.
        pt = {0:0.15, 4:0.2}
    n_frames: int
        Number of frames required

    Returns
    -------
       list of tuples where first number is the trajectory index and
       second is the frame index
    """
    dimensions = list(pt_dict.keys())
    d_data = [i[:, dimensions] for i in data]

    tree = KDTree(d_data)
    pt = [pt_dict[i] for i in dimensions]
    dis, ind = tree.query(pt, n_frames)
    return ind

def find_nearest(a, a0, prev_pt=None):

    b= a[:,:,tic_index]
    "Element in nd array `a` closest to the scalar value `a0`"
    dis_to_wanted = np.abs(b - a0)
    if prev_pt is None:
        idx = np.nanargmin(dis_to_wanted)
    else:
        #need to min distance along
        sorted_dis = np.argsort(dis_to_wanted)

    return b.flat[idx]

def _frame_loader(yaml_file, prt, key_list, indices, save_trj, fname=None):
    traj_list=[]
    #if not an iterator make it
    if type(indices)!=list and indices.shape==(2,):
        indices = [indices]
    for ind in indices:
        traj_index, frame_index = ind
        traj_name = key_list[traj_index]
        traj_list.append(load_frame(yaml_file["base_dir"],
                                    prt.name, yaml_file["protein_dir"],
                                    traj_name, frame_index))

    trj = traj_list[0] + traj_list[1:]
    if save_trj and fname is not None:
        with enter_protein_mdl_dir(yaml_file, prt.name):
            trj.save_xtc("%s"%fname)
            if not os.path.isfile("prot.pdb"):
                trj[0].save_pdb("prot.pdb")

    return trj

def pull_frames(yaml_file, prt, tic_index, n_frames,scheme="linear", save_trj=True):
    """
    :param yaml_file: The loaded yaml file
    :param prt The prtein mdl.
    :param tic_index: tic index to sample along
    :param n_frames:number of watned frames
    :param scheme:One of 3
    linear:Samples the tic linearly
    random:Samples the tic randomly
    edge: Samples the tic edges only
    :return:
    :output: This will write out a log file and a xtc file. The log file will
    contain the values of the tic that were obtained while the xtc file will contain
    the tic itself. Also returns the trajector obj
    """

    #get some statistics about the data
    key_list = list(prt.tica_data.keys())
    data = [prt.tica_data[i] for i in key_list]
    #get indices
    indices = sample_dimension(data, tic_index, n_frames, scheme)
    #load traj
    fname = "tic%d.xtc"%tic_index
    trj =_frame_loader(yaml_file, prt, key_list, indices, save_trj,fname)

    return trj

def _load_protein_matrices(yaml_file, protein_name):
    """
    Helper routine to load matrices for a protein
    :param yaml_file: yaml file to work with
    :param protein_name: name of the protein
    :return:
     prj :The protein Series
     prt : The protein project
     key_mapping: mapping of the assigment matrix 0-axis to traj names
     assignment_matrix: Massive matrix of
     tics_mapping: mapping of the tics_array matrix 0-axis to traj names
     tics_array: Multi dimensional array where the 0th axis is equal to the
     number of trajectors, the 1st axis is equal to largest traj and the
     3rd dimension is equal to the number of tics in the mdl.
    """
    prj = ProteinSeries(yaml_file)
    prt = Protein(prj, protein_name)

    key_mapping, assignment_matrix  = create_assignment_matrix(prt.fixed_assignments)
    tics_mapping , tics_array  = create_tics_array(prt.fixed_assignments, prt.kmeans_mdl,
                                        prt.tica_data)

    return prj, prt, key_mapping, assignment_matrix, tics_mapping, tics_array

def sample_one_tic(yaml_file,protein_name,tic_index,n_frames, scheme="linear"):
    """
    :param yaml_file: The project's yaml file
    :param protein: The name of protein
    :param tic_index: Tic index to sample along
    :param n_frames: The number of frames wanted
    :return: Dumps a tic%d.xtc and tic%d.log for a given
    protein inside its model.
    """
    yaml_file = load_yaml_file(yaml_file)
    prj = ProteinSeries(yaml_file)
    prt = Protein(prj, protein_name)

    return pull_frames(yaml_file, prt, tic_index, n_frames, scheme)

def sample_tic_region(yaml_file, protein_name, tic_region,
                      n_frames=50, fname=None,save_trj=True):
    """
    Helper function for sampling tic in a particular tic_region.
    :param yaml_file: The projects yaml file
    :param protein_name: The name of the protein
    :param tic_region(dict): The tic_region. Can be multidimensional with
    1 number per tic coordinate(defaults to 0 for all non-mentioned regions)
    :param n_frames: The number of frames around the coordinate
    :return:
    """

    yaml_file = load_yaml_file(yaml_file)

    prj = ProteinSeries(yaml_file)
    prt = Protein(prj, protein_name)

    key_list = list(prt.tica_data.keys())
    data = [prt.tica_data[i] for i in key_list]
    indices = sample_region(data, tic_region, n_frames)

    if fname is None:
        fname = "sampled_tic_region.xtc"
    trj =_frame_loader(yaml_file, prt, key_list, indices, save_trj, fname)

    return trj


def sample_for_all_proteins(yaml_file, protein=None, tics=None, n_frames=100,
                            scheme="linear"):
    """
    :param yaml_file: The project yaml file.
    :param protein: The name of the protein. If none, then it is
    done for all the protein names in the yaml_file. If it is a list,
    it is iteratively done for each of the protein else its only called
    once.
    :param tics: list of tics to sample from. If None, then
    it is done for all the tics specified in the yaml file
    :param n_frames number of frames wanted for each tic
    :param scheme:One of 3 sampling schemes
    linear:Samples the tic linearly
    random:Samples the tic randomly
    edge: Samples the tic edges only
    :return:
    """

    yaml_file = load_yaml_file(yaml_file)
    if protein is None :
        protein =  yaml_file[protein_list]

    if tics==None:
        tics = range(yaml_file["params"]["tica__n_components"])

    for protein_name in protein:
        for tic_index in tics:
            sample_one_tic(yaml_file, protein_name, tic_index, n_frames,
                           scheme)

    return


def _map_tic_component(tic_component, df, trj):
    '''
    Function map a tic component to all atoms and optionally all residues
    by summing over all the feature importances where the residue index
    appears
    :param tic_component: The feature weight vector to use
    :param df: Dataframe describing each feature
    :param trj: mdtraj trajctory obj
    :return:atom_importance and residue importance
    '''
    n_atoms = trj.top.n_atoms
    n_residues = trj.top.n_residues

    atom_importance_vector = np.zeros((1,n_atoms))
    residue_importance_vector =  np.zeros((1,n_residues))

    #over all features
    for i in df.iterrows():
        #over all every residue in that feature
        for j in i[1]["resids"]:
            #add to running total
            residue_importance_vector[0, j] += abs(tic_component[i[0]])

    for r in trj.topology.residues:
        for a in r.atoms:
            atom_importance_vector[0, a.index] = residue_importance_vector[0, r.index]


    return atom_importance_vector, residue_importance_vector

def max_movement(tica_data, index=0,num_wanted=1):
    """
    Helper routine to find trajectory that goes the max along any coordinate
    :param tica_data: Tica data
    :param index: tIC index , defaults to 0
    :param number of trajectories wanted
    :return: list of trajectory names where the first traj
    has the highest movement in that dimension
    """
    diff_list=[]
    key_list=[]
    for i in tica_data.keys():
        diff_list.append(np.max(tica_data[i][:,index]) - np.min(tica_data[i][:,index]))
        key_list.append(i)

    return [key_list[i] for i in np.argsort(diff_list)[-num_wanted:]][::-1]
