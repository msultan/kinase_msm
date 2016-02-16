from kinase_msm.feature_selection import series_feature_slicer
from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import enter_protein_data_dir
from kinase_msm.series_setup import setup_series_analysis
from test_pipeline import create_fake_data
import os, glob
from msmbuilder.utils import verboseload, verbosedump


def test_slicer():
    with enter_temp_directory():
        base_dir = os.path.abspath(os.path.curdir)
        mdl_dir = os.path.join(base_dir,"mdl_dir")
        feature_dir = "feature_dir"
        series_name = "fake_series"
        protein_list = ["kinase_1", "kinase_2"]
        project_dict = {"kinase_1": ["fake_proj1",],
                        "kinase_2": ["fake_proj2"]}
        mdl_params = {'tica__n_components': 1, 'tica__lag_time': 1,
                  'tica__weighted_transform': True, 'tica__shrinkage': 0.01,
                  'cluster__n_clusters': 2,
                  'msm__lag_time': 1, 'bootstrap__n_samples':1 }

        create_fake_data(base_dir, protein_list, project_dict)

        yaml_file = setup_series_analysis(base_dir, mdl_dir, feature_dir,
                                  series_name, protein_list,
                                  project_dict, mdl_params)

        dict_feat_ind={}
        dict_feat_ind["kinase_1"] =[0, 2]
        dict_feat_ind["kinase_2"] =[1, 1, 0, 2]

        series_feature_slicer(yaml_file, dict_feat_ind)


        for protein in protein_list:
            with enter_protein_data_dir(yaml_file, protein):
                assert (os.path.isdir("sliced_feature_dir"))
                flist = glob.glob("./%s/*.jl"%feature_dir)
                for fname in flist:
                    original_file = verboseload(fname)
                    expected_file = original_file[:, dict_feat_ind[protein]]
                    written_file = verboseload("./%s/%s"%("sliced_feature_dir",
                                                          os.path.basename(fname)
                                                          ))
                    assert (expected_file==written_file).all()
    return