#!/bin/env/python
__author__ = 'muneeb'
import numpy as np


def create_assignment_matrix(assignments):
    """
    :param assignments: Dictionary of assignments keyed by the trajectory name.
    :return: key_mapping: Dictionary that maps the array index to assignments file name
    :return: assignment_matrix: matrix of assignments of shape n_traj, length of max traj.
    Non-used array positions are set to -1 . Useful for doing fast assignment lookups
    """
    n_traj = np.shape(assignments.keys())[0]
    max_length = np.max([np.shape(assignments[i]) for i in assignments.keys()])
    assignment_matrix = np.zeros((n_traj, max_length)) - 1
    key_mapping = {}
    for i, v in enumerate(assignments.keys()):
        current_traj_length = np.shape(assignments[v])[0]
        assignment_matrix[i, :current_traj_length] = assignments[v]
        key_mapping[i] = v

    return key_mapping, assignment_matrix


def create_tics_array(assignments, kmeans_mdl, tica_data):
    """
    :param assignments: Dictionary of assignments keyed by the trajectory name.
    :param kmeans_mdl: The KMeans mdl used to build the state decomposition.
    :param tica_data: The tica data keyed by trajectory name
    :return: tics_array: 3-D matrix os shape n_traj, max_length, num_of_tics.
    Useful for doing fast tic lookups. All non-used values are NANs
    """
    n_traj = np.shape(assignments.keys())[0]
    max_length = np.max([np.shape(assignments[i]) for i in assignments.keys()])

    tics_to_use = kmeans_mdl.cluster_centers_.shape[1]
    tics_array = np.empty((n_traj, max_length, tics_to_use))
    tics_array[:, :, :] = np.NAN
    for index, trjname in enumerate(assignments.keys()):
        current_traj_length = np.shape(assignments[trjname])[0]
        for tic_index in range(tics_to_use):
            tics_array[index, :current_traj_length, tic_index] = tica_data[trjname][:, tic_index]
    return tics_array