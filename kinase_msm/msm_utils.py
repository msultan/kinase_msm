#!/bin/evn python
import numpy as np
import os
from msmbuilder.utils import verbosedump
from multiprocessing import cpu_count, Pool
from .data_loader import load_yaml_file
from .data_loader import load_frame, enter_protein_mdl_dir
from .data_transformer import create_assignment_matrix
from .mdl_analysis import ProteinSeries, Protein
from msmbuilder.utils.nearest import KDTree
from .tica_utils import _frame_loader

def _sample_state(jt):
    #get where the state exists
    state, assignment_matrix, key_mapping, base_dir, prt_name = jt
    trj, frame = np.where(assignment_matrix==state)
    index = np.random.randint(0,len(trj))
    filename = key_mapping[trj[index]]
    frame_index = frame[index]
    trj_frame = load_frame(base_dir, prt_name, filename, frame_index)
    return trj_frame


def sample_discarded_states(yaml_file, prt_list=None):
    """
    :param yaml_file: The model yaml file to work with
    :param prt_list:
    :return:
    """
    raise NotImplementedError("Sorry :(")

def sample_states(yaml_file, prt_name):
    yaml_file = load_yaml_file(yaml_file)
    ser = ProteinSeries(yaml_file)
    prt =Protein(ser, prt_name)

    return

def kdtree_maker(prt):
    key_list = list(prt.tica_data.keys())
    data = [prt.tica_data[i] for i in key_list]

    tree = KDTree(data)
    return key_list, tree

def sample_state_centroid(yaml_file, prt_name, states='all', n_frames=10, output_name="centroids.xtc"):
    """
    :param yaml_file: The series' yaml file
    :param prt_name: The name of the protein
    :param state: Int or list of ints or 'all'
    :param n_frames: Number of output frames
    :return: output name(FILE IS DUMPED AS XTC)
    """
    yaml_file = load_yaml_file(yaml_file)
    ser = ProteinSeries(yaml_file)
    prt =Protein(ser, prt_name)


    if type(states)==str and states=='all':
        states = np.arange(prt.n_states_)
    elif type(states)==list or type(states)==int:
        states = np.array(states)
    else:
        raise ValueError("States needs to  either a int or list or all")

    key_list, tree = kdtree_maker(prt)

    trj_list=[]
    for state in states:
        inv_mapped_state = prt.msm.inverse_transform([state])[0][0]
        kcent = prt.kmeans_mdl.cluster_centers_[inv_mapped_state,:]
        dis, ind = tree.query(kcent, k=n_frames)
        trj_list.append(_frame_loader(yaml_file, prt, key_list,
                                      indices=ind, save_trj=False, fname=None))

    trj = trj_list[0] +trj_list[1:]
    with enter_protein_mdl_dir(yaml_file, prt_name):
        trj.save_xtc(output_name)

    return trj


def sample_msm_traj(yaml_file, prt_name, n_steps, starting_state = None,
                    fname="msm_traj.xtc"):
    """
    :param yaml_file: The model's yaml file
    :param prt: The name of the protein mdl
    :param n_steps: The number of markovian frames desired.
    :param starting_state: If None, we start from the most populated state.
    :param fname : The output filename
    :return: Dumps the msm traj.
    """

    yaml_file = load_yaml_file(yaml_file)
    ser = ProteinSeries(yaml_file)
    prt =Protein(ser,prt_name)

    # this returns in original assignment space
    msm_traj = prt.msm.sample_discrete(state=starting_state, n_steps=n_steps)
    # there we use the original assignment matrix too
    key_mapping, assignment_matrix = create_assignment_matrix(prt.assignments)

    jbs =[(state, assignment_matrix, key_mapping, ser.base_dir, prt.name) for state in msm_traj]
    p = Pool(int(cpu_count()/4))
    trj_list = p.map(_sample_state, jbs)
    print("Done")
    trj = trj_list[0] + trj_list[1:]

    with enter_protein_mdl_dir(yaml_file, prt_name):
        verbosedump(msm_traj,"msm_traj.pkl")
        trj.save_xtc(fname)
        if not os.path.isfile("prot.pdb"):
            trj[0].save_pdb("prot.pdb")
    return
