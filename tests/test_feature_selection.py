from kinase_msm.feature_selection import series_feature_slicer, _map_residue_ind_seq_ind
from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import enter_protein_data_dir, load_yaml_file, load_random_traj
from kinase_msm.series_setup import setup_series_analysis
from test_pipeline import create_fake_data
from nose.tools import with_setup
from test_convert_series import _setup_test, _cleanup_test
import os, glob
from msmbuilder.utils import verboseload, verbosedump

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))


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

@with_setup(_setup_test, _cleanup_test)
def test_map_residue_seq_no_insert():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)

        expected[protein] = [i.index for i in t.top.residues]
        aligned_dict[protein] = [i.code for i in t.top.residues]

        actual=_map_residue_ind_seq_ind(yaml_file, protein, aligned_dict)

        assert expected[protein] == list(actual.values())

    return

def test_map_residue_seq_with_insert():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)

        expected[protein] = [i.index+3 for i in t.top.residues]
        aligned_dict[protein] = list("___")+ [i.code for i in t.top.residues]

        actual=_map_residue_ind_seq_ind(yaml_file, protein, aligned_dict)
        assert expected[protein] == list(actual.values())

    return

def test_map_residue_seq_with_insert_after_10():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)
        #add an insertion AFTER 10 residues. We expect all but the 10 have

        expected[protein] = [i for i in range(10)] + [i+3 for i in range(10, t.n_residues)]

        aligned_dict[protein] = [i.code for i in t.top.residues][:10]+\
                                list("___")+ \
                                [i.code for i in t.top.residues][10:]

        actual=_map_residue_ind_seq_ind(yaml_file, protein, aligned_dict)
        assert expected[protein] == list(actual.values())

    return

def test_map_residue_seq_with_insert_at_end():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)
        #add an insertion AFTER 10 residues. We expect all but the 10 have

        expected[protein] = [i for i in range(t.n_residues)]

        aligned_dict[protein] = [i.code for i in t.top.residues]+list("___")

        actual=_map_residue_ind_seq_ind(yaml_file, protein, aligned_dict)

        assert expected[protein] == list(actual.values())
    return


def test_map_residue_seq_with_two_inserts():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)
        #add an insertion AFTER 10 residues. and then again at 20
        expected[protein] = [i for i in range(10)] + [i+3 for i in range(10, 20)]+\
                            [i+5 for i in range(20, t.n_residues)]

        aligned_dict[protein] = [i.code for i in t.top.residues][:10]+\
                                list("___")+ \
                                [i.code for i in t.top.residues][10:20]+\
                                list("__")+ \
                                [i.code for i in t.top.residues][20:]

        actual=_map_residue_ind_seq_ind(yaml_file, protein, aligned_dict)
        assert expected[protein] == list(actual.values())

    return


