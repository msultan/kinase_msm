#!/bin/evn python

from kinase_msm.tica_utils import *
import pandas as pd
import mdtraj as mdt
from kinase_msm.fit_transform_kinase_series import *
from kinase_msm.msm_utils import sample_msm_traj
from kinase_msm.tica_utils import _map_tic_component
from kinase_msm.mdl_analysis import ProteinSeries, Protein
from test_convert_series import _setup_test, _cleanup_test
from nose.tools import with_setup
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.utils import verboseload


if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))


def _test_tic_sampling(yaml_file, protein_name, tic_list, n_frames, scheme):
    #test to make sure we are sampling right
    sample_for_all_proteins(yaml_file, [protein_name],
                            tic_list, n_frames, scheme=scheme)
    ser = ProteinSeries(yaml_file)
    prt = Protein(ser, protein_name)

    for tic_index in [0,1]:
        traj_path = os.path.join(base_dir,yaml_file["mdl_dir"],
                                 protein_name,"tic%d.xtc"%tic_index)
        traj_top = os.path.join(base_dir,yaml_file["mdl_dir"],
                                protein_name, "prot.pdb")
        tica_traj = mdt.load(traj_path,top=traj_top)
        feat = DihedralFeaturizer(types=['phi', 'psi','chi1'])

        f = feat.partial_transform(tica_traj)
        t_f = np.round(prt.tica_mdl.transform([f]))

        #check that the tic goes from min to max
        assert t_f[0][0][tic_index] <= t_f[0][-1][tic_index]
        all_vals = []
        for traj_tica_data in prt.tica_data.values():
            all_vals.extend(traj_tica_data[:,tic_index])
            #sort it because all three sampling schemes use it

        all_vals = np.round(np.sort(all_vals))
        print(tic_index)
        print(t_f[0][:,tic_index] >= all_vals[0])
        print(t_f[0][:,tic_index] <= all_vals[-1])
        #make sure the frames are within limitsss
        assert (t_f[0][:,tic_index] >= all_vals[0]).all()
        assert (t_f[0][:,tic_index] <= all_vals[-1]).all()
    return True


def _test_sample_region(yaml_file, protein_name, tic_region,
                      n_frames=5, fname="temp.xtc"):
    sample_tic_region(yaml_file, protein_name, tic_region,
                      n_frames, fname)
    assert os.path.isfile(os.path.join(yaml_file["mdl_dir"],
                                              protein_name,fname))
    return True

def test_tica_utils():
    np.random.seed(42)
    yaml_file = os.path.join(base_dir,"mdl_dir","project.yaml")
    yaml_file = load_yaml_file(yaml_file)
    fit_pipeline(yaml_file["base_dir"])
    assert _test_tic_sampling(yaml_file, "kinase_1", [0,1], 5, "linear")
    a={}
    a[0]=0.3
    a[1]=0.4
    assert _test_sample_region(yaml_file,"kinase_1", a)
    return True

def _fit_transform(prt, trj):
    f=DihedralFeaturizer(types=['phi', 'psi','chi1'])
    feat = f.partial_transform(trj)
    t_f = prt.tica_mdl.transform([feat])
    st = prt.kmeans_mdl.transform(t_f)
    return st

def test_msm_traj():
    yaml_file = os.path.join(base_dir,"mdl_dir","project.yaml")
    yaml_file = load_yaml_file(yaml_file)
    n_steps=2
    ser = ProteinSeries(yaml_file,base_dir)
    prt = Protein(ser, "kinase_2")
    starting_state = prt.msm.state_labels_[0]
    sample_msm_traj(yaml_file, "kinase_2",n_steps=n_steps,starting_state=starting_state)
    with enter_protein_mdl_dir(yaml_file, "kinase_2"):
        msm_steps = verboseload("msm_traj.pkl")
        msm_traj = mdt.load("msm_traj.xtc",top="prot.pdb")
        assert (msm_traj.n_frames==n_steps)
        assert(len(msm_steps)==n_steps)
        states = _fit_transform(prt, msm_traj)
        assert (states==msm_steps).all()


@with_setup(_setup_test, _cleanup_test)
def test_map_tic_component():
    yaml_file = os.path.join(base_dir,"mdl_dir","project.yaml")
    yaml_file = load_yaml_file(yaml_file)
    fit_pipeline(yaml_file["base_dir"])

    with enter_protein_data_dir(yaml_file, "kinase_1"):
        df = pd.DataFrame(verboseload(
            os.path.join(yaml_file["feature_dir"],
                         "feature_descriptor.h5")
        ))
        trj = mdt.load(os.path.join("protein_traj", "fake_proj1_0_0.hdf5"))


    ser = ProteinSeries(yaml_file,base_dir)
    prt = Protein(ser, "kinase_1")

    tica_mdl = prt.tica_mdl
    tic_index=0
    t_c = tica_mdl.components_[tic_index, :]

    a_i, r_i = _map_tic_component(t_c, df, trj)

    assert len(a_i[0]) == trj.n_atoms
    assert len(r_i[0]) == trj.n_residues

    #spot check residue 0
    df2 = pd.DataFrame([i[1] for i in df.iterrows() if (i[1]["resid"]==0).any()])
    r0_imp = np.sum(abs(t_c[df2.index]))
    assert r0_imp==r_i[0,0]
