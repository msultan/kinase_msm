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
series_name : {series_name}
kinase_list : {kinase_list}
project_dict : {project_dict}
mdl_params : {mdl_params}
"""


def setup_series_analysis(base_dir, mdl_dir, feature_dir, series_name, kinase_list,
                          project_dict, mdl_params):
    """
    :param base_dir: Directory where all the data is found
    :param series_name: Series name for the set of kinases
    :param kinase_list: List of kinases
    :param project_dict: dictionary of lists where each list
    holds the projects corresponding to the kinase
    :param mdl_params: Mdl that is being built
    :return:The script will backup and current mdl_dir and make a new
    directory with a clean set of subfolders and a yaml file. The
    yaml file will also be returned to directly be fed into the fit transform.
    """
    if not os.path.isdir(base_dir):
        raise Exception("Base directory doesn't exist")

    # make sure all the data exists
    for kinase in kinase_list:
        if not os.path.isdir(os.path.join(base_dir, kinase)):
            raise Exception("Directory for kinase %s doesn't exist" % kinase)
        if kinase not in project_dict.keys():
            raise Exception("Projects for kinase %s doesn't exist" % kinase)
        for project in project_dict[kinase]:
            if not os.path.isdir(os.path.join(base_dir, kinase, project)):
                raise Exception("Project %s for kinase %s doesn't exist"
                                % (project, kinase))

    if os.path.isdir(mdl_dir):
        ctime = str(int(time.time()))
        shutil.move(mdl_dir, mdl_dir + ctime)
        warnings.warn("Moved Previous mdl dir to %s%s" % (mdl_dir,ctime))
    else:
        pass

    os.mkdir(mdl_dir)
    for kinase in kinase_list:
        os.mkdir(os.path.join(mdl_dir, kinase))

    with open(os.path.join(mdl_dir, 'project.yaml'), 'w') as yaml_out:
        yaml_file = yaml.load(yaml_template.format(base_dir=base_dir,
                                                                 mdl_dir=mdl_dir,
                                                                 feature_dir=feature_dir,
                                                                 series_name=series_name,
                                                                 kinase_list=kinase_list,
                                                                 project_dict=project_dict,
                                                                 mdl_params=mdl_params))
        yaml_out.write(yaml.dump(yaml_file))

    return yaml_file
