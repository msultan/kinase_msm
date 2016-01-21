#!/bin/evn python

from kinase_msm.tica_utils import *
import mdtraj as mdt
from kinase_msm.fit_transform_kinase_series import *

from msmbuilder.featurizer import DihedralFeaturizer

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
        t_f = prt.tica_mdl.transform([f])

        #check that the tic goes from min to max
        assert t_f[0][0][tic_index] <= t_f[0][-1][tic_index]
        all_vals = []
        for traj_tica_data in prt.tica_data.values():
            all_vals.extend(traj_tica_data[:,tic_index])
            #sort it because all three sampling schemes use it

        all_vals = np.sort(all_vals)
        print(tic_index)
        print(t_f[0][:,tic_index] >= all_vals[0])
        print(t_f[0][:,tic_index] <= all_vals[-1])
        #make sure the frames are within limitsss
        assert (t_f[0][:,tic_index] >= all_vals[0]).all()
        assert (t_f[0][:,tic_index] <= all_vals[-1]).all()
    return False


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
    assert _test_sample_region(yaml_file,"kinase_1",a)
    return True

