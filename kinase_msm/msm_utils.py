#!/bin/evn python
from kinase_msm.data_loader import load_yaml_file
import numpy as np
import os
from msmbuilder.utils import verbosedump
from multiprocessing import Pool, cpu_count
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from kinase_msm.data_loader import load_frame, load_current_protein_model, enter_protein_mdl_dir
from kinase_msm.data_transformer import create_assignment_matrix


def _sample_state(jt):
    #get where the state exists
    state, assignment_matrix, key_mapping, base_dir, prt_name = jt
    trj,frame = np.where(assignment_matrix==state)
    while True:
        try:
            index = np.random.randint(0,len(trj))
            filename = key_mapping[trj[index]]
            frame_index = frame[index]
            frame = load_frame(base_dir, prt, filename, frame_index)
            break
        except:
            continue
    return frame

    return

def sample_discarded_states(yaml_file, prt_list=None):
    """
    :param yaml_file: The model yaml file to work with
    :param prt_list:
    :return:
    """
    raise NotImplementedError("Sorry :(")

def sample_states(yaml_file, prt):

    return


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


    base_dir, mdl_dir, msm_mdl, \
    tica_mdl, tica_data, kmeans_mdl, \
    assignments = load_current_protein_model(yaml_file, prt_name)

    with enter_protein_mdl_dir(yaml_file, prt_name):

        if starting_state is None:
            starting_state = np.argmax(msm_mdl.populations_)

        msm_traj = msm_mdl.sample_discrete(state=starting_state,n_steps=n_steps)
        verbosedump(msm_traj,"msm_traj.pkl")

        key_mapping, assignment_matrix = create_assignment_matrix(assignments)

        p=Pool(cpu_count()/4)
        jbs =[(state, assignment_matrix, key_mapping, base_dir, prt_name) for state in msm_traj]
        trj_list=p.map(_sample_state, jbs)

        trj = trj_list[0] + trj_list[1:]
        trj.save_xtc(fname)
        if not os.path.isfile("prot.pdb"):
            trj[0].save_pdb("prot.pdb")


    return
