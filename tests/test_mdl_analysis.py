#!/bin/env python
import os
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from kinase_msm.fit_transform_kinase_series import fit_pipeline
from kinase_msm.series_setup import setup_series_analysis
from msmbuilder.decomposition import tICA
from msmbuilder.msm import MarkovStateModel
from msmbuilder.utils import verboseload
from mdtraj.utils.contextmanagers import enter_temp_directory
from test_pipeline import create_fake_data
import numpy as np

def test_project():
    with enter_temp_directory():
        base_dir = os.path.abspath(os.path.curdir)
        print(base_dir)
        print(type(base_dir))
        mdl_dir = os.path.join(base_dir,"mdl_dir")
        feature_dir = "feature_dir"
        series_name = "fake_series"
        protein_list = ["kinase_1", "kinase_2"]
        project_dict = {"kinase_1": ["fake_proj1",],
                        "kinase_2": ["fake_proj2"]}
        mdl_params = {'tica__n_components': 1, 'tica__lag_time': 1,
                  'tica__weighted_transform': True, 'tica__shrinkage': 0.01,
                  'cluster__n_clusters': 2,'msm__lag_time': 1, 'bayesmsm__n_samples':1,
                  'bayesmsm__n_steps':1}

        create_fake_data(base_dir, protein_list, project_dict)

        setup_series_analysis(base_dir, mdl_dir, feature_dir,
                                  series_name, protein_list,
                                  project_dict, mdl_params)
        fit_pipeline(base_dir)

        prj = ProteinSeries(os.path.join(mdl_dir,"project.yaml"))

        assert isinstance(prj, ProteinSeries)
        assert isinstance(prj.tica_mdl ,tICA)

        assert _test_protein_without_project()
        assert _test_protein_with_project(prj)
        assert _test_tic_dict(prj)

    return


def _test_protein_without_project():
    p1 = None
    try:
        p1 = Protein("fake_series", "kinase_1")
    except:
        assert not isinstance(p1, Protein)
        return True
    return

def _test_tic_dict(prj):
    p1 = Protein(prj, "kinase_1")
    #p1._get_all_tics()
    current_data = []
    rnd_tic = np.random.randint(p1.n_tics_)
    rnd_state = np.random.randint(p1.n_states_)
    for traj_index, traj_name in enumerate(p1.fixed_assignments.keys()):
        for f_i, f in enumerate(p1.fixed_assignments[traj_name]):
            if f == rnd_state:
                current_data.append(p1.tica_data[traj_name][f_i][rnd_tic])
    assert len(current_data)==len(p1.tic_data(rnd_tic)[rnd_state])
    assert current_data == p1.tic_data(rnd_tic)[rnd_state]
    return True

def _test_protein_with_project(prj):
    p1 = Protein(prj, "kinase_1")
    p2 = Protein(prj, "kinase_2")
    assert isinstance(p1, Protein)
    assert isinstance(p1.msm, MarkovStateModel)
    assert (p1.msm.left_eigenvectors_ ==
            verboseload(os.path.join(prj.mdl_dir,"kinase_1","msm_mdl.pkl")).left_eigenvectors_).all()

    assert (p2.msm.left_eigenvectors_ ==
            verboseload(os.path.join(prj.mdl_dir,"kinase_2","msm_mdl.pkl")).left_eigenvectors_).all()

    return True