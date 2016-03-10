from kinase_msm.feature_selection import series_feature_slicer, \
    _map_residue_ind_seq_ind,_present_for_all, _get_common_residues, _get_common_features
from mdtraj.utils.contextmanagers import enter_temp_directory
from kinase_msm.data_loader import enter_protein_data_dir, load_yaml_file, load_random_traj
from kinase_msm.series_setup import setup_series_analysis
from test_pipeline import create_fake_data
import os, glob
from msmbuilder.utils import verboseload, verbosedump
from msmbuilder.featurizer import DihedralFeaturizer

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

        expected[protein] = [i.index for i in t.top.residues if i.is_protein]
        aligned_dict[protein] = t.top.to_fasta(chain=0)
        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
        print(list(actual.values()))
        print(expected[protein])
        assert expected[protein] == list(actual.values())

    return

def test_map_residue_seq_with_insert():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        expected = {}

        t = load_random_traj(yaml_file, protein)

        expected[protein] = [i.index+3 for i in t.top.residues if i.is_protein]
        aligned_dict[protein] = "---"+ t.top.to_fasta(chain=0)
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

        expected[protein] = [i for i in range(10) if t.top.residue(i).code is not None] + \
                            [i+3 for i in range(10, t.n_residues) if t.top.residue(i).code is not None]
        prt_code = t.top.to_fasta(chain=0)
        aligned_dict[protein] = prt_code[:10]+\
                                "---"+ \
                                prt_code[10:]
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

        expected[protein] = [i for i in range(t.n_residues) if t.top.residue(i).code is not None]

        aligned_dict[protein] = t.top.to_fasta(chain=0)+"---"

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
        expected[protein] = [i for i in range(10) if t.top.residue(i).code is not None] + \
                            [i+3 for i in range(10, 20)  if t.top.residue(i).code is not None]+\
                            [i+5 for i in range(20, t.n_residues) if t.top.residue(i).code is not None]

        prt_code = t.top.to_fasta(chain=0)
        aligned_dict[protein] = prt_code[:10]+\
                                "---"+ \
                                prt_code[10:20]+\
                                "--"+ \
                                prt_code[20:]

        aligned_seq = aligned_dict[protein]
        actual,_ =_map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
        assert expected[protein] == list(actual.values())

    return

def test_present_for_all_same_seq():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        aligned_dict[protein] = t.top.to_fasta(chain=0)

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

    aligned_dict["p1"] = ["A","S","D","B","-","-","A","S","D"]
    aligned_dict["p2"] = ["A","M","D","B","-","M","A","S","D"]

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

    aligned_dict["p1"] = ["-","A","S","D","B","-","-","A","S","D"]
    aligned_dict["p2"] = ["-","A","M","D","B","-","M","A","S","D"]

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

    aligned_dict["p1"] = ["-","A","S","D","B","-","-","A","S","X"]
    aligned_dict["p2"] = ["-","A","M","D","B","-","M","A","S","D"]

    prt_mapping["p1"] ={0:1, 1:2, 2:3, 3:4, 4:7, 5:8,6:9}
    prt_mapping["p2"] ={0:1, 1:2, 2:3, 3:4, 4:6, 5:7,6:8,7:9}

    p1_res=_present_for_all("p1", prt_mapping["p1"], prt_seq["p1"], aligned_dict)
    p2_res=_present_for_all("p2", prt_mapping["p2"], prt_seq["p2"], aligned_dict)

    assert(p1_res==[0,2,3,4,5])
    assert(p2_res==[0,2,3,5,6])
    return

def test_get_common_residues():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        aligned_dict[protein] = t.top.to_fasta(chain=0)

    res_dic =  _get_common_residues(yaml_file, aligned_dict)
    for protein in yaml_file["protein_list"]:
        print(len(res_dic[protein]),t.n_residues)
        assert(len(res_dic[protein])==len(t.top.to_fasta(chain=0)))

    return

def test_get_common_features():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        aligned_dict[protein] = t.top.to_fasta(chain=0)

    common_res_dic =  _get_common_residues(yaml_file, aligned_dict)

    f= DihedralFeaturizer()
    common_feature_dic = _get_common_features(yaml_file,f, common_res_dic, False)
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        assert(len(common_feature_dic[protein])==f.transform(t)[0].shape[1])

    return


def test_get_common_features_2():
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    aligned_dict={}
    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        aligned_dict[protein] = t.top.to_fasta(chain=0)

    common_res_dic =  _get_common_residues(yaml_file, aligned_dict)
    #lets get rid of last few residues
    for protein in yaml_file["protein_list"]:
        common_res_dic[protein] = common_res_dic[protein][:5]
        print(common_res_dic[protein])

    f= DihedralFeaturizer(types=['phi','psi','chi1'])
    common_feature_dic = _get_common_features(yaml_file,f, common_res_dic, False)

    for protein in yaml_file["protein_list"]:
        t = load_random_traj(yaml_file, protein)
        # 4 phi, 4 psi and 4 chi1(one of the 5 has no chi1).
        # times 2 to account for sincos transoform
        assert(len(common_feature_dic[protein])==(4+4+4)*2)

    return