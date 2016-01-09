#!/bin/evn python
from kinase_msm.data_loader import load_yaml_file
import numpy as np
from kinase_msm.mdl_analysis import Project, Protein
from kinase_msm.data_loader import load_frame
from kinase_msm.data_transformer import create_assignment_matrix, create_tics_array

"""
Set of helper scripts for sampling tica
"""

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.nanargmin(np.abs(a - a0))
    return a.flat[idx]


def pull_frames(yaml_file, protein_name, tic_index, n_frames, key_mapping,
                     assignment_matrix, tics_array,tica_data,scheme="linear"):
    """
    :param yaml_file: The loaded yaml file
    :param protein_name: name of the protein
    :param tic_index: tic index to sample along
    :param n_frames:number of watned frames
    :param key_mapping:mapping of len of matrix of assignments to the
     traj names
    :param assignment_matrix:matrix of assignment
    :param tics_array:3d array of all tica daata
    :param tica_data:Dictionary of tica files
    :param scheme:Only linearly sampled
    :return:
    :output: This will write out a log file and a xtc file. The log file will
    contain the values of the tic that were obtained while the xtc file will contain
    the tic itself.
    """

    #get some statistics about the data
    max_tic_movement = max([max(i[:,tic_index]) for i in tica_data.values()])
    min_tic_movement = min([min(i[:,tic_index]) for i in tica_data.values()])

    #get lineraly placed points
    if scheme=="linear":
        lin_place_points = np.linspace(min_tic_movement, max_tic_movement, n_frames)
    traj_list = []
    actual_tic_val_list=[]
    for v,i in enumerate(lin_place_points):
        actual_tic_val = find_nearest(tics_array[:,:,tic_index],i)
        actual_tic_val_list.append(actual_tic_val)

        traj_index, frame_index = np.where(tics_array[:,:,tic_index]==actual_tic_val)
        traj_name = key_mapping[traj_index[0]]
        actual_tic_val_list.append([i, actual_tic_val,traj_name,frame_index[0]])
        traj_list.append(load_frame(yaml_file["base_dir"],
                                    protein_name,traj_name,frame_index[0]))

    trj = traj_list[0]
    for i in traj_list[1:]:
        trj += i

    save_dir = os.path.join(yaml_file["mdl_dir"],protein_name)
    #dump the log file
    with open(os.path.join("tic%d.log"%tic_index,"w")) as fout:
        fout.write("Tic Value, Actual Value, TrajName, FrmInd\n")
        for line in actual_tic_val_list:
            fout.write("%s\n"%line)
    trj.save_xtc(os.path.join(save_dir,"tic%d.xtc"%tic_index))

    trj[0].save_pdb(os.path.join(save_dir,"prot.pdb"))

    return


def sample_one_tic(yaml_file,protein_name,tic_index,n_frames, scheme="linear"):
    """
    :param yaml_file: The project's yaml file
    :param protein: The name of protein
    :param tic_index: Tic index to sample along
    :param n_frames: The number of frames wanted
    :return: Dumps a tic%d.xtc and tic%d.log for a given
    protein inside its model.
    """
    prj = Project(yaml_file)
    prt = Protein(prj, protein_name)

    key_mapping, assignment_matrix  = create_assignment_matrix(prt.fixed_assignments)
    _ , tics_array  = create_tics_array(prt.tica_data)

    yaml_file = load_yaml_file(yaml_file)
    pull_frames(yaml_file,protein_name,tic_index,n_frames,key_mapping,assignment_matrix,
                tics_array, prt.tica_data,scheme)
    return


def sample_all_tics(yaml_file, protein_name, n_frames, scheme="linear"):
    """
    :param yaml_file: The project's yaml file
    :param protein: The name of protein
    :param n_frames: The number of frames needed for each tic
    :return: Dumps the tic%d_lin.xtc in the protein mdl
    folders for each of the tics
    """
    prj = Project(yaml_file)
    prt = Protein(prj, protein_name)

    for tic_index in range(prt.n_tics_):
        sample_one_tic(yaml_file, protein_name, tic_index, n_frames, scheme)

    return


def sample_for_all_proteins(yaml_file, protein=None, tics=[0], n_frames=100,
                            scheme="linear"):
    """
    :param yaml_file: The project yaml file.
    :param protein: The name of the protein. If none, then it is
    done for all the protein names in the yaml_file. If it is a list,
    it is iteratively done for each of the protein else its only called
    once.
    :param tics: list of tics to sample from. If None, then
    it is only done for the dominant tic
    :return:
    """
    yaml_file = load_yaml_file(yaml_file)
    if protein is None :
        for protein in yaml_file[protein_list]:
            sample_all_tics(yaml_file, protein, n_frames)
    else:
        for protein_name in protein:
            sample_all_tics(yaml_file, protein_name, n_frames)



    return
