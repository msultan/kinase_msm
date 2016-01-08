#!/bin/evn python
import mdtraj as mdt
from kinase_msm.data_loader import load_yaml_file
import numpy as np
from data_loader import change_protein_data_dir, change_protein_mdl_dir
"""
Set of helper scripts for doing tica analysis

"""

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.nanargmin(np.abs(a - a0))
    return a.flat[idx]


def pull_tica_frames(mutant,tic_index,n_frames,key_mapping,assignment_matrix,
                     tics_array,tica_data,save_dir):
    change_mutant_dir(mutant)

    #a few nice things to pre calculate
    n_traj = assignment_matrix.shape[0]
    max_traj_len =assignment_matrix.shape[1]

    #figure out how much trajectories move through tic space in general

    #get some statistics about the data
    max_tic_movement = max([max(i[:,tic_index]) for i in tica_data.values()])
    min_tic_movement = min([min(i[:,tic_index]) for i in tica_data.values()])

    lin_place_points = np.linspace(min_tic_movement,max_tic_movement,n_frames)
    traj_list = []
    actual_tic_val_list=[]
    for v,i in enumerate(lin_place_points):
        actual_tic_val = find_nearest(tics_array[:,:,tic_index],i)
        actual_tic_val_list.append(actual_tic_val)

        traj_index,frame_index = np.where(tics_array[:,:,tic_index]==actual_tic_val)
        traj_name = key_mapping[traj_index[0]]

        #print v,i, actual_tic_val, traj_name,frame_index[0]
        traj_list.append(load_frame(mutant,traj_name,frame_index[0]))

    trj = traj_list[0]
    for i in traj_list[1:]:
        trj += i

    trj.save_xtc(os.path.join(save_dir,"tic%d_lin_sampled.xtc"%tic_index))
    trj[0].save_pdb(os.path.join(save_dir,"prot.pdb"))
    return


def sample_along_tic(yaml_file, protein=None, tics=[0], n_frames=100):
    """
    :param yaml_file: The project yaml file.
    :param protein: The name of the protein. If none, then it is
    done for all the protein names in the yaml_file
    :param tics: list of tics to sample from. If None, then
    it is only done for the dominant tic
    :return:
    """
    yaml_file = load_yaml_file(os.path.join(mdl_dir, "project.yaml"))
    key_mapping,assignment_matrix = create_assignment_matrix(assignments)
    tics_array = create_tics_array(assignments, kmeans_mdl, tica_data)

    base_dir,wt_dir,msm_mdl,tica_mdl,tica_data,kmeans_mdl,assignments = \
        load_current_protein_model(mutant,sanity=True)

    #for tic_index in tics:


    return
