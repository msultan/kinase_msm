# kinase_msm
[![Build Status](https://travis-ci.com/msultan/kinase_msm.svg?token=Qs64kEbR7UUaswHXtepX&branch=master)](https://travis-ci.com/msultan/kinase_msm)

set of scripts to analyze kinases

These scripts are designed to work with the following directory structure. The base directory has all the
kinases associated with that set of projects. The kinase project
base_dir
project.yaml
|->kinase1
   |-> kinase1_proj_1
       |-> trajectories
   |-> kinase1_proj_2
       |-> trajectories
   |-> allprojects_features_dir
   |-> protein_traj
|->kinase2
   |-> kinase2_proj1
       |->trajectories
   |-> allprojects_features_dir
   |-> protein_traj
It will generate a mdl_dir with sub_dir for each kinase along with other data in the mdl dir.
|->mdl_dir
   |-> kinase1
   |-> kinase2

