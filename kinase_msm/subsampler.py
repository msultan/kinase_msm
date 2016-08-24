from kinase_msm.data_loader import load_yaml_file, enter_protein_data_dir
from kinase_msm.featurize_project import _check_output_folder_exists
import yaml
import os,glob
import mdtraj as md
import itertools
from multiprocessing import Pool, cpu_count

def subsample_traj(jt):
    inp_file,output_file,stride = jt
    t = md.load(inp_file, stride=stride)
    t.save_hdf5(output_file)
    return

def subsample_protein(yaml_file, protein, stride=5,out_dir="sub_protein_traj"):
    yaml_file=load_yaml_file(yaml_file)

    p=Pool(int(cpu_count()/2))

    with enter_protein_data_dir(yaml_file, protein):
        flist = [os.path.abspath(i) for i in
                 glob.glob("%s/*.hdf5"%yaml_file["protein_dir"])]

    base_dir = yaml_file["base_dir"]
    new_output_dir = os.path.join(base_dir,protein,out_dir)
    if not os.path.isdir(new_output_dir):
        os.mkdir(new_output_dir)
    fout = [os.path.join(new_output_dir,os.path.basename(i)) for i in flist]

    zippy = zip(flist, fout, itertools.repeat(stride))

    jobs= [(i,o,s) for i,o,s in zippy]
    p.map(subsample_traj,jobs)
    return

def subsample_series(yaml_file,stride=5,out_dir="sub_protein_traj",overwrite=True):
    yaml_file = load_yaml_file(yaml_file)
    for protein in yaml_file["protein_list"]:
        subsample_protein(yaml_file,protein, stride, out_dir)
    yaml_file["protein_dir"] = out_dir
    #write the new yaml file
    if overwrite:
        with open(os.path.join(yaml_file["mdl_dir"],
                           'project.yaml'), 'w') as yaml_out:
            yaml_out.write(yaml.dump(yaml_file))

    return
