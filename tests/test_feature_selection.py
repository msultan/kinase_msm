from kinase_msm.feature_selection import series_feature_slicer, _map_residue_ind_seq_ind,_present_for_all
from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import enter_protein_data_dir, load_yaml_file, load_random_traj
from kinase_msm.series_setup import setup_series_analysis
from test_pipeline import create_fake_data
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


def test_map_residue_seq_no_insert():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)

        expected[protein] = [i.index for i in t.top.residues]
        aligned_dict[protein] = [i.code for i in t.top.residues]
        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)

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
        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
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
        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
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

        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)

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

        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
        assert expected[protein] == list(actual.values())

    return

def test_present_for_all_same_seq():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        aligned_dict[protein] = [i.code for i in t.top.residues]

    for protein in yaml_file["protein_list"]:
        aligned_seq = aligned_dict[protein]
        prt_mapping, prt_seq =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
        assert(len(_present_for_all(protein, prt_mapping, prt_seq, aligned_dict))==len(prt_seq))
    return


def test_present_for_all_2():

    aligned_dict={}
    prt_seq ={}
    prt_mapping={}

    prt_seq["p1"] = ["A","S","D","B","A","S","D"]
    prt_seq["p2"] = ["A","M","D","B","M","A","S","D"]

    aligned_dict["p1"] = ["A","S","D","B","_","_","A","S","D"]
    aligned_dict["p2"] = ["A","M","D","B","_","M","A","S","D"]

    prt_mapping["p1"] ={0:0,1:1, 2:2, 3:3, 4:6, 5:7,6:8}
    prt_mapping["p2"] ={0:0,1:1, 2:2, 3:3, 4:5, 5:6,6:7, 7:8}

    p1_res=_present_for_all("p1", prt_mapping["p1"], prt_seq["p1"], aligned_dict)
    p2_res=_present_for_all("p2", prt_mapping["p2"], prt_seq["p2"], aligned_dict)

    print(p1_res)
    assert(p1_res==[0,2,3,4,5,6])
    assert(p2_res==[0,2,3,5,6,7])
    return


def test_present_for_all_3():

    aligned_dict={}
    prt_seq ={}
    prt_mapping={}

    prt_seq["p1"] = ["A","S","D","B","A","S","D"]
    prt_seq["p2"] = ["A","M","D","B","M","A","S","D"]

    aligned_dict["p1"] = ["_","A","S","D","B","_","_","A","S","D"]
    aligned_dict["p2"] = ["_","A","M","D","B","_","M","A","S","D"]

    prt_mapping["p1"] ={0:1,1:2, 2:3, 3:4, 4:7, 5:8,6:9}
    prt_mapping["p2"] ={0:1,1:2, 2:3, 3:4, 4:6, 5:7,6:8,7:9}

    p1_res=_present_for_all("p1", prt_mapping["p1"], prt_seq["p1"], aligned_dict)
    p2_res=_present_for_all("p2", prt_mapping["p2"], prt_seq["p2"], aligned_dict)

    assert(p1_res==[0,2,3,4,5,6])
    assert(p2_res==[0,2,3,5,6,7])
    return



def test_present_for_all_4():

    aligned_dict={}
    prt_seq ={}
    prt_mapping={}

    prt_seq["p1"] = ["A","S","D","B","A","S","X"]
    prt_seq["p2"] = ["A","M","D","B","M","A","S","D"]

    aligned_dict["p1"] = ["_","A","S","D","B","_","_","A","S","X"]
    aligned_dict["p2"] = ["_","A","M","D","B","_","M","A","S","D"]

    prt_mapping["p1"] ={0:1, 1:2, 2:3, 3:4, 4:7, 5:8,6:9}
    prt_mapping["p2"] ={0:1, 1:2, 2:3, 3:4, 4:6, 5:7,6:8,7:9}

    p1_res=_present_for_all("p1", prt_mapping["p1"], prt_seq["p1"], aligned_dict)
    p2_res=_present_for_all("p2", prt_mapping["p2"], prt_seq["p2"], aligned_dict)

    assert(p1_res==[0,2,3,4,5])
    assert(p2_res==[0,2,3,5,6])
    return


