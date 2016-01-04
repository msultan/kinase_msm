#!/bin/env/python
import numpy as np
from kinase_msm.data_transformer import create_assignment_matrix, \
    create_tics_array


def test_create_assignment_matrix():
    assignments = {}
    assignments["traj1"] = [0, 1, 2]
    assignments["traj2"] = [2, 1, 2]
    assignments["traj3"] = [2, 1, 0, 3]

    key_mapping, assignment_matrix = create_assignment_matrix(assignments)

    assert assignment_matrix.shape == (3, 4)
    assert (assignment_matrix[1][0:3] == assignments[key_mapping[1]][:3]).all()

    return


class fake_kmeans():
    def __init__(self, n_states, n_tics):
        self.cluster_centers_ = np.zeros((n_states, n_tics))


def test_create_tics_array():
    n_traj = 100
    max_traj_length = 100
    n_states = 10
    n_tics = 5

    test_kmeans_mdl = fake_kmeans(n_states, n_tics)
    tica_data = {}
    test_assignments = {}

    for i in range(n_traj):
        test_assignments[str(i)] = np.random.randint(n_states,
                                                     size=max_traj_length)
        tica_data[str(i)] = np.random.normal(size=(max_traj_length,
                                                   n_tics))

    key_mapping, tics_array = create_tics_array(test_assignments,
                                                test_kmeans_mdl,
                                                tica_data)

    assert tics_array.shape == (n_traj, max_traj_length, n_tics)
    for i in range(10):
        temp_num = np.random.randint(n_traj)
        assert (tics_array[temp_num] == tica_data[key_mapping[temp_num]]).all()
    return
