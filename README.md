# kinase_msm
[![Build Status](https://travis-ci.com/msultan/kinase_msm.svg?token=Qs64kEbR7UUaswHXtepX&branch=master)](https://travis-ci.com/msultan/kinase_msm)

set of scripts to analyze kinases

These scripts are designed to work with the following directory structure. The base directory has all the
kinases associated with that set of projects. The kinase project \n
+ base_dir 
+ project.yaml 
+ kinase1
  + kinase1_proj_1 
     + trajectories 
  + kinase1_proj_2
     + trajectories 
  + features_dir 
  + protein_only_traj 
+ kinase2
  + kinase2_proj_1 
     + trajectories 
  + features_dir 
  + protein_only_traj 

It will generate a mdl_dir in the base_dir with sub directories for each kinase along with other data in the mdl dir. The general idea is to compare multiple kinases with each other or to make model for a single kinase only.  
+ mdl_dir 
   + kinase1
   + kinase2
