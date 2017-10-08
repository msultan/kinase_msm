#!/bin/env python
import mdtraj as mdt
from kinase_msm.data_loader import  load_yaml_file
from kinase_msm.subsampler import subsample_series
import os,glob

if os.path.isdir("tests"):
    base_dir = os.path.abspath(os.path.join("./tests/test_data"))
else:
    base_dir = os.path.abspath(os.path.join("./test_data"))

def test_subsampler():
    print(base_dir)
    yaml_file = load_yaml_file(os.path.join(base_dir,"mdl_dir","project.yaml"))
    out_dir = "sub_protein_traj"
    subsample_series(yaml_file,out_dir=out_dir,overwrite=False)
    assert(os.path.isdir(os.path.join(base_dir,"kinase_1",out_dir)))
    for k in ["kinase_1","kinase_2"]:
        for i in glob.glob(os.path.join(base_dir,k, "protein_traj","*.hdf5")):
            t1 = mdt.load(i)
            t2 = mdt.load(os.path.join(base_dir,k ,out_dir, os.path.basename(i)))
            assert (t1.n_frames==t2.n_frames*5)