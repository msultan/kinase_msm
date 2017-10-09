#!/bin/env python

from __future__ import print_function
import os
import yaml
import warnings
import shutil
import time

yaml_template = """
base_dir : {base_dir}
mdl_dir : {mdl_dir}
feature_dir: {feature_dir}
protein_dir: {protein_dir}
series_name : {series_name}
protein_list : {protein_list}
project_dict : {project_dict}
mdl_params : {mdl_params}
"""


def setup_series_analysis(base_dir, mdl_dir, feature_dir, series_name, protein_list,
                          project_dict,mdl_params=None, protein_dir="protein_traj"):
    """
    This script sets up a framework for directories to perform multiple analyses of a series of
    related proteins.  This may be a database of mutants, a family of related proteins, etc.

	Input:

    :param base_dir: Directory where all the data is found. base_dir should be organized as

        base_dir
          > protein_dir
            > project_dir

    where protein_dir corresponds to the names of proteins in the series, and project_dir
    corresponds to data for a specific analysis.  TODO: don't understand project delineation here?

    :param series_name: Series name for the series of proteins.

    :param protein_list: List of proteins in the series.

    :param project_dict: Dictionary keyed by the names of the proteins in the series, 
          where each entry is a list of projects corresponding to that protein.

    :param mdl_params: Parameters for the model that is being fit (tICA, MSM, etc). 

    :return: The script will create the directory mdl_dir with the following structure:

         mdl_dir
           > protein_dir
             > project_dir
           project.yaml

     where project.yaml details the model being fit.  If a directory named mdl_dir currently exists,
    it will create a backup. 
    TODO: ? The yaml file will also be returned to directly be fed into the fit transform.
    """
    if not os.path.isdir(base_dir):
        raise Exception("Base directory doesn't exist")

    # make sure all the data exists
    for protein in protein_list:
        if not os.path.isdir(os.path.join(base_dir, protein)):
            raise Exception("Directory for protein %s doesn't exist" % protein)
        if protein not in project_dict.keys():
            raise Exception("Projects for protein %s doesn't exist" % protein)
        for project in project_dict[protein]:
            if not os.path.isdir(os.path.join(base_dir, protein, project)):
                raise Exception("Project %s for protein %s doesn't exist"
                                % (project, protein))

# Write yaml file in base_dir, detailing series analysis setup.

    with open(os.path.join(base_dir,"series.yaml"),'w') as yaml_out:
        yaml_file = yaml.load(yaml_template.format(base_dir=base_dir,
                                                   mdl_dir=mdl_dir,
                                                   feature_dir=feature_dir,
                                                   series_name=series_name,
                                                   protein_dir=protein_dir,
                                                   protein_list=protein_list,
                                                   project_dict=project_dict,
                                                   mdl_params=None
                                                   ))
        yaml_out.write(yaml.dump(yaml_file))

    if mdl_params is not None:
        if os.path.isdir(mdl_dir):
            ctime = str(int(time.time()))
            shutil.move(mdl_dir, mdl_dir + ctime)
            warnings.warn("Moved Previous mdl dir to %s%s" % (mdl_dir, ctime))
        else:
            pass

# make mdl_dir that has subdirs for each protein in the series.
        os.mkdir(mdl_dir)
        for protein in protein_list:
            os.mkdir(os.path.join(mdl_dir, protein))

# write project.yaml file in mdl_dir, detailing project setup.
        with open(os.path.join(mdl_dir, 'project.yaml'), 'w') as yaml_out:
            yaml_file["mdl_params"]=mdl_params
            yaml_out.write(yaml.dump(yaml_file))

    return yaml_file
