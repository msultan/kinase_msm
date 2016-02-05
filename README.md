# kinase_msm
[![Build Status](https://travis-ci.com/msultan/kinase_msm.svg?token=Qs64kEbR7UUaswHXtepX&branch=master)](https://travis-ci.com/msultan/kinase_msm)

set of scripts to analyze kinases

These scripts are designed to work with the following directory structure. The base directory has all the
kinases associated with that set of projects. The kinase project \n
+ base_dir 
+ series.yaml 
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
   + project.yaml 
   + kinase1
   + kinase2


Usecase:
``` python
from kinase_msm.series_setup import setup_series_analysis
from kinase_msm.fit_transform_kinase_series import fit_pipeline
import os 

base_dir = "/hsgs/nobackup/msultan/research/kinase/btk_kinase/fah_data/rcsb"

protein_list =["rcsb_mdl", "rcsb_hmdl"]

proj_dict={}

proj_dict={}
proj_dict["rcsb_mdl"] =["PROJ9145","PROJ9151"]
proj_dict["rcsb_hmdl"] = ["PROJ9146","PROJ9152"]

#using her2 mdl osprey results
osprey_stride = 100 

osprey_mdl_params = {'tica__n_components': 5, 'tica__lag_time': 3, 'tica__weighted_transform': True, 'tica__shrinkage': None, 'cluster__n_clusters': 500}


mdl_params = osprey_mdl_params
mdl_params["tica__lag_time"] = mdl_params["tica__lag_time"]*osprey_stride
mdl_params["msm__lag_time"] = 3*osprey_stride

mdl_dir = os.path.join(base_dir,"mdl_dir")
feature_dir = "feature_dir"
series_name = "btk_series"
project_dict = proj_dict

def main():
    yaml_file = setup_series_analysis(base_dir,mdl_dir,feature_dir,series_name,protein_list, proj_dict, mdl_params)
    featurize_series(yaml_file)
    fit_pipeline(base_dir)

if __name__=="__main__":
    main()

```
