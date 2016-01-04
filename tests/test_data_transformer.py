#!/bin/env/python
from kinase_msm.data_transformer import create_assignment_matrix, create_tics_array


def test_create_assignment_matrix():
    test_assignments = {}
    test_assignments["traj1"] = [0, 1, 2]
    test_assignments["traj2"] = [2, 1, 2]
    test_assignments["traj3"] = [2, 1, 0, 3]

    key_mapping, assignment_matrix = create_assignment_matrix(test_assignments)

    assert key_mapping[0] == "traj1"
    assert key_mapping[2] == "traj3"
    assert assignment_matrix.shape == (3,4)
    assert (assignment_matrix[1][0:3] == test_assignments["traj2"][:3]).all()

    return

def test_create_tics_array():
    return
