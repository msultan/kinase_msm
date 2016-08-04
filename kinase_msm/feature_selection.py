from Bio import SeqIO
import glob
import os
import numpy as np
import pandas as pd
from multiprocessing import Pool
from msmbuilder.utils import verboseload, verbosedump
from kinase_msm.data_loader import load_yaml_file
from kinase_msm.featurize_project import _check_output_folder_exists
from kinase_msm.data_loader import load_random_traj, \
    enter_protein_data_dir, enter_protein_mdl_dir
from sklearn.base import clone
from msmbuilder.featurizer import ContactFeaturizer
import itertools
"""
Set of routines to select common features amongst
proteins based upon their sequence similarity.
Creates a new folder to dump features of interest in there.
"""

def _parse_alignment_file(filename):
    aligned_dict = {}
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    for fasta in fasta_sequences:
        aligned_dict[fasta.id]=fasta.seq.tostring().upper()
    return aligned_dict


def _present_for_all(protein, prt_mapping, prt_seq, aligned_dict):
    #test to make sure the seq matches up
    aligned_seq = aligned_dict[protein]
    prt_seq = ''.join(prt_seq)
    assert(prt_seq==''.join([i for i in aligned_seq if i!="-"]))

    result_vector = []
    #for every index and code

    max_length = max([len(i) for i in aligned_dict.values()])

    #for that position in the alignment
    position = 0
    while position < max_length:
        possible_codes = set([aligned_dict[p][position] for p in aligned_dict.keys()])
        print(possible_codes)
        # if we find a insertion, we skip the that and the next residue
        if "-" in possible_codes:
            position += 1
            continue
        # if only one code there.
        prev_codes = set([aligned_dict[p][position-1] for p in aligned_dict.keys()])
        if len(possible_codes)==1 and ("-" not in prev_codes):
            index = [key for key in prt_mapping.keys() if prt_mapping[key]==position][0]
            result_vector.append(index)
        position += 1
    return result_vector


def _map_residue_ind_seq_ind(yaml_file, protein, aligned_seq, trj=None):
    if trj is None:
        trj = load_random_traj(yaml_file, protein)
    mapping = {}
    seq_index = 0
    prt_seq = trj.top.to_fasta(chain=0)
    #test to make sure the alignment sequence matches with the protein sequence.
    #get rid of _ from the alignment to account for additions/deletions.
    assert(prt_seq==''.join([i for i in aligned_seq if i!="-"]))
    for i in [i.index for i in trj.top.residues if i.is_protein]:
        while True:
            if trj.top.residue(i).code == aligned_seq[seq_index]:
                mapping[i] = seq_index
                seq_index += 1
                break
            else:
                seq_index += 1
                continue

    return mapping, prt_seq


def _get_common_residues(yaml_file, aligned_dict):
    #for every protein in the list
    result_dict={}
    for protein in yaml_file["protein_list"]:
        #get the mapping
        aligned_seq = aligned_dict[protein]
        prt_mapping, prt_seq = _map_residue_ind_seq_ind(yaml_file, protein, aligned_seq)
        result_dict[protein] = _present_for_all(protein, prt_mapping, prt_seq, aligned_dict)
    return result_dict, prt_mapping


def _get_common_features(yaml_file, featurizer, aligned_dict,save_df=True):
    """
    Function to get the common features across protein using the common residues.
    can optionally save the pandas data to the mdl_dir
    :param yaml_file: The protein yaml_file
    :param featurizer: featurizer object used.
    :param prt_mapping: Mapping of each residue to its sequence
    :param aligned_dict : Dictionary of alignments for each protein
    :return:
    """
    result_dict = {}
    df_dict={}
    for protein in yaml_file["protein_list"]:
        print(protein)
        #reset the featurizer
        featurizer = clone(featurizer)
        trj = load_random_traj(yaml_file, protein)
        df = pd.DataFrame(featurizer.describe_features(trj))
        prt_mapping, prt_seq = _map_residue_ind_seq_ind(yaml_file, protein,
                                                        aligned_dict[protein])
        feature_vec =[]
        #for every feature
        for i in df.iterrows():
            #get the index and the feature itself
            feature_ind, feature_dict = i
            all_res_in_algn = []
            mapped_index_list=[]
            for aa_ind in feature_dict["resids"]:
                aa_code = prt_seq[aa_ind]
                #make sure we have the same residue
                assert(trj.top.residue(aa_ind).code==aa_code)
                #get the mapping for that aa to the main alignment
                mapped_index = prt_mapping[aa_ind]
                #for every protein in the alignment, check if we have the same residue
                #at the same position
                all_res_in_algn.append(np.alltrue([aligned_dict[prt][mapped_index]==aa_code
                                          for prt in yaml_file["protein_list"]]))
                mapped_index_list.append(mapped_index)


            #to account for additions and deletions, we check if the difference between
            #the mapping and the actual residue codes is the same.
            mapped_index_difference = [x - mapped_index_list[i - 1]
                                       for i, x in enumerate(mapped_index_list) if i > 0]
            resid_index_difference = [x - feature_dict["resids"][i - 1]
                                       for i, x in enumerate(feature_dict["resids"]) if i > 0]
            if not np.all(mapped_index_difference==resid_index_difference):
                all_res_in_algn.append(False)


            if np.alltrue(all_res_in_algn):
                feature_vec.append(feature_ind)

        df_dict[protein] = df.iloc[feature_vec]
        result_dict[protein] = feature_vec

        if save_df:
            new_df = df.iloc[feature_vec]
            with enter_protein_mdl_dir(yaml_file, protein):
                verbosedump(new_df, os.path.join("feature_descriptor.h5"))
            with enter_protein_data_dir(yaml_file, protein):
                verbosedump(new_df, os.path.join("sliced_feature_dir",
                                                 "feature_descriptor.h5"))
    return result_dict, df_dict

def _slice_file(job_tuple):
    inp_file, feature_ind, output_folder = job_tuple
    featurized_file = verboseload(inp_file)
    sliced_file = featurized_file[:, feature_ind]
    sliced_file_out = os.path.join(output_folder, os.path.basename(inp_file))
    verbosedump(sliced_file, sliced_file_out)
    return

def _feature_slicer(yaml_file, dict_feat_ind, folder_name, view):
    protein_list = yaml_file["protein_list"]


    for protein in protein_list:
        _check_output_folder_exists(yaml_file, protein, folder_name)
        feature_folder = os.path.join(yaml_file["base_dir"],
                                  protein, yaml_file["feature_dir"])
        output_folder = os.path.join(yaml_file["base_dir"],
                                  protein, folder_name)
        flist = glob.glob(os.path.join(feature_folder,"*.jl"))

        feature_ind = dict_feat_ind[protein]
        jobs = [(inp_file, feature_ind, output_folder) for inp_file in flist]
        view.map(_slice_file, jobs)

    return

def series_feature_slicer(yaml_file, dict_feat_ind=None,
                          featurizer=None,
                          folder_name="sliced_feature_dir",
                         view=None):

    """
    :param yaml_file: The project yaml file with
    :param dict_feat_ind: Dict of wanted feature indices for each protein. Defaults to
    none when you want the code to figure out what features to keep.
    :param featurizer: The featurizer object that was used to generat.
    :param folder_name: Name of the output folder. Defaults to sliced_feature_dir
    :param view: pool of workers. Defaults to multiprocessing
    :return: None
    """

    yaml_file = load_yaml_file(yaml_file)

    if view is None:
        view = Pool()

    #if we want to do this and we cant find the sequence
    if dict_feat_ind is None and ("alignment_file" not in yaml_file
                                  or featurizer is None
                                  or (not hasattr(featurizer, "describe_features"))):
        raise ValueError("To find common features, we need both "
                         "the alignment file in the yaml file"
                         "AND a featurizer obj that supports describe_features")


    if dict_feat_ind is None:
        #load alignment file
        aligned_dict = _parse_alignment_file(yaml_file["alignment_file"])
        #get list of common residue indices
        #dict_common_res, prt_mapping = _get_common_residues(yaml_file, aligned_dict)
        #get list of feature indices
        dict_feat_ind, df_dict = _get_common_features(yaml_file, featurizer, aligned_dict)

    _feature_slicer(yaml_file, dict_feat_ind, folder_name, view)

    return



def test_series_slicer(yaml_file, folder_name="sliced_feature_dir"):
    yaml_file = load_yaml_file(yaml_file)

    df_dict={}
    for protein in yaml_file["protein_list"]:
        with enter_protein_data_dir(yaml_file, protein):
            df_dict[protein] = verboseload(os.path.join(os.getcwd(),
                folder_name,"feature_descriptor.h5"))
    for ind,protein in enumerate(yaml_file["protein_list"]):
        for ind2, protein2 in enumerate(yaml_file["protein_list"]):
            assert (df_dict[protein].resnames==
                    df_dict[protein2].resnames).all()

    return


def create_equivalent_contact_featurizer(yaml_file, alignment_file,
                                         protein_list=None,
                                         wanted_sequence_ind_locs=None,
                                         same_residue=True,
                                         **kwargs):
    """
    Create a equivalent contacts featurizer for a set of proteins
    :param yaml_file: yaml file location
    :param alignment_file: alignment file location
    :param wanted_sequence_ind_locs: wanted sequence index positions in the alignment
    You need to just figure out the wanted location for one residue.
    _map_residue_ind_seq_ind function can help with this
    :same residue: True is you would restrict to having the same residue at the same
    sequence position.
    :param kwargs: kwargs for the contact featurizer
    :return: dictionary of contact featurizers. one for each protein
    """
    featurizer_dict={}

    #load alignment file
    yaml_file = load_yaml_file(yaml_file)
    alignment_file = _parse_alignment_file(alignment_file)
    if protein_list is None:
        protein_list = yaml_file["protein_list"]

    if wanted_sequence_ind_locs is None:
        #use the max length(probably a horrible idea)
        max_seq_len = max([len(alignment_file[i]) for i in alignment_file.keys()])
        wanted_sequence_ind_locs = [i for i in range(max_seq_len)]

    for protein in protein_list:
        print(protein)
        #get a list of residues we can keep
        can_keep=[]
        #get mapping and seq
        prt_mapping, prt_seq = _map_residue_ind_seq_ind(yaml_file, protein,
                                                        alignment_file[protein])
        #for wanted positions in the massive wanted indices list
        inv_map = {v: k for k, v in prt_mapping.items()}

        for position in wanted_sequence_ind_locs:
            #get the
            #get the possible codes at every position
            possible_codes = set([alignment_file[p][position] for p in alignment_file.keys()])
            #if there is not a missing residue

            if not "-" in possible_codes:
                if same_residue and len(set(possible_codes))!=1:
                    continue
                # get the inverse mapping and add it to the list of can keep
                residue_index = inv_map[position]
                can_keep.append(residue_index)
        #sort it because i dont want random bs issues.
        can_keep = np.sort(can_keep)
        #get its pairs
        pairs = [i for i in itertools.combinations(can_keep, 2)]

        featurizer_dict[protein] = ContactFeaturizer(contacts=pairs, **kwargs)

    return featurizer_dict
