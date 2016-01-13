#!/bin/evn python

import os
from msmbuilder.utils import verboseload, verbosedump
from kinase_msm.data_loader import load_yaml_file


class Project(object):

    def __init__(self, yaml_file):
        self.yaml_file = load_yaml_file(yaml_file)
        self.base_dir = self.yaml_file["base_dir"]
        self.mdl_dir = self.yaml_file["mdl_dir"]
        self.relative_loc = os.path.split(yaml_file)[0]
        self.kmeans_mdl = verboseload(
            os.path.join(self.relative_loc, "kmeans_mdl.pkl"))
        self.tica_mdl = verboseload(os.path.join(self.relative_loc, "tica_mdl.pkl"))


class Protein(object):
    """
    protein class to load all msm related data for a single protein within a project
    """

    def __init__(self, series, name):
        if not isinstance(series, Project):
            raise Exception("We need a project series to be associated "
                            "with this kinase")
        self.name = name
        self.project = series
        self.kmeans_mdl = self.project.kmeans_mdl
        self.tica_mdl = self.project.tica_mdl
        self.protein_mdl_dir = os.path.join(self.project.relative_loc, self.name)
        self.msm = verboseload("%s/msm_mdl.pkl" % self.protein_mdl_dir)
        self.bayes_mdl = verboseload(
            '%s/bayesmsm_mdl.pkl' % self.protein_mdl_dir)
        self.tica_data = verboseload("%s/tica_data.pkl" % self.protein_mdl_dir)
        self.assignments = verboseload(
            "%s/assignments.pkl" % self.protein_mdl_dir)
        self.fixed_assignments = verboseload(
            "%s/fixed_assignments.pkl" % self.protein_mdl_dir)
        self.n_states_ = self.msm.n_states_
        self.n_tics_ = self.kmeans_mdl.cluster_centers_.shape[1]
        self.tic_dict = {}

    def _get_all_tics(self):
        """
        fill in all the tics at once in the tic dict which is useful
        :return:
        """
        if not self.tic_dict:
            for tic_index in range(self.n_tics_):
                self.tic_dict[tic_index] = {}
                for j in range(self.n_states_):
                    self.tic_dict[tic_index][j] = []

            for traj_index, traj_name in enumerate(self.fixed_assignments.keys()):
                # for all the fixed_state assignments
                for f_i, fixed_state in enumerate(self.fixed_assignments[traj_name]):
                    # go through and find the
                    try:
                        for tic_index in range(self.n_tics_):
                            self.tic_dict[tic_index][fixed_state].append(
                                self.tica_data[traj_name][f_i][tic_index])
                    except:
                        pass
            return
        else:
            return

    def tic_data(self, tic_index):
        """
        Helper function to get a dictionary keyed on the state with all the
        tic values for that state in it.
        :param tic_index:
        :return: a dictionary of lists
        """
        self._get_all_tics()
        return self.tic_dict[tic_index]
