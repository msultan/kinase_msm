#!/bin/evn python

import pandas as pd
import numpy as np

"""
set of helper routines to plot things
"""

def global_tic_boundaries(prt_list, tic_list, n_bins=100):

    """
    :param prt_list: list of proteins loaded using the Protein class
    :param tic_list: list of tics to compute for
    :return:a dictionary where the key is the tic index and the value is
    list containing the linearly spaced tic value going from the global min
    and max
    """
    assert isinstance(prt_list,list)
    results_dict = {}
    for tic_index in tic_list:
        min_tic_val = min([prot.tic_min[tic_index] for prot in prt_list])
        max_tic_val =  max([prot.tic_max[tic_index] for prot in prt_list])
        results_dict[tic_index] = np.linspace(min_tic_val,max_tic_val,n_bins)

    return results_dict

def histogram_data(prj, prt, tic_list, x_array=None, y_array=None, n_bins=100):

    #simple check for making sure we have a list
    if not isinstance(tic_list, list):
        tic_list = [tic_list]

    c_x = prt.tic_dict[tic_list[0]]

    if x_array is None:
        lin_spaced_tic = global_tic_boundaries([prt],tic_list,n_bins)
        x_array = lin_spaced_tic[tic_list[0]]
        if y_array is None and len(tic_list)==2:
            y_array =lin_spaced_tic[tic_list[1]]
    else:
        n_bins = len(x_array)-1

    x_center = (x_array[:-1] + x_array[1:]) / 2

    if len(tic_list)==1:
        H_overall = np.zeros(n_bins)

    #if y is also needed
    if len(tic_list)==2:
        c_y = prt.tic_dict[tic_list[1]]
        H_overall = np.zeros((n_bins,n_bins))

    H={}
    for i in range(prt.n_states_):
        if len(tic_list)==1:
            H[i],x=np.histogram(c_x[i], bins=x_array, normed=True)
        elif len(tic_list)==2:
            H[i],x,y=np.histogram2d(c_x[i], c_y[i], bins=[x_array, y_array], normed=True)
        else:
            raise Exception("cant do this")

        H_overall = H_overall + prt.msm.populations_[i]*H[i]

    #convert to kcal/mol by
    H_overall = -0.6*np.log(H_overall)

    return H, H_overall.transpose(), x_center

def bayes_one_dim_free_energy(prj,prt,tic_index,n_bins=100 ,lin_spaced_tic=None):
    """
    :param prj: Project that the protein is a part of
    :param prt: the protein itself
    :param tic_index: The tic index that is needed
    :param n_bins: n_bins
    :param lin_spaced_tic: linearly sampled tic
    :param errorbars: whether or not to compute multiple free enegies
    using bayes msm mdl
    :return: a pandas dataframe containing the free energy for the
    every point along the tic coordinate. The mdl index column contains
    an integer for each of the bayesian mdls.
    """
    free_energy = []

    if lin_spaced_tic is None:
        lin_spaced_tic = global_tic_boundaries([prt],tic_index,n_bins)[tic_index]
    else:
        n_bins = len(lin_spaced_tic) - 1

    #get the centers stacked nicely
    tic_cen = np.repeat([(lin_spaced_tic[:-1] + lin_spaced_tic[1:]) / 2],
                        prt.bayes_mdl.n_samples,axis=0).flatten()
    protein_name = np.repeat(prt.name,prt.bayes_mdl.n_samples*n_bins).flatten()

    mdl_index = np.array([np.repeat(i,n_bins)
                          for i in range(prt.bayes_mdl.n_samples)]).flatten()

    #get data
    H,H_msm,_ = histogram_data(prj,prt,tic_index,x_array=lin_spaced_tic)

    for i in range(prt.bayes_mdl.n_samples):
        H_overall=np.zeros(n_bins)
        for j in range(prt.n_states_):
            H_overall = H_overall + prt.bayes_mdl.all_populations_[i][j]*H[j]
        #convert to free enenrgy
        H_overall = -0.6*np.log(H_overall)

        free_energy.extend(H_overall)

    df=pd.DataFrame([list(tic_cen),list(free_energy),list(protein_name),list(mdl_index)]).T
    df.columns=["tic_value","free_energy","protein_name","mdl_index"]
    return df

def one_dim_free_energy(prj, prt, tic_index, n_bins=100 ,
                        lin_spaced_tic=None, errorbars=False):
    """
    :param prj: Project that the protein is a part of
    :param prt: the protein itself
    :param tic_index: The tic index that is needed
    :param n_bins: n_bins
    :param lin_spaced_tic: linearly sampled tic
    :param errorbars: whether or not to compute multiple free enegies
    using bayes msm mdl
    :return: a pandas dataframe containing the free energy for the
    every point along the tic coordinate. The mdl index column contains
    "mle" for the msm free energy and an integer for the bayesian mdl
    """
    free_energy = []

    if lin_spaced_tic is None:
        lin_spaced_tic = global_tic_boundaries([prt],tic_index,n_bins)[tic_index]
    else:
        n_bins = len(lin_spaced_tic) - 1

    #get the centers stacked nicely
    tic_center = np.repeat([(lin_spaced_tic[:-1] + lin_spaced_tic[1:]) / 2],
                           1 , axis=0).flatten()
    protein_name = np.repeat(prt.name,n_bins).flatten()

    mdl_index = np.repeat("mle",n_bins).flatten()

    #get data
    H,H_msm,_ = histogram_data(prj,prt,[tic_index],x_array=lin_spaced_tic)

    free_energy.extend(H_msm)

    msm_df=pd.DataFrame([list(tic_center),list(free_energy),
                         list(protein_name),list(mdl_index)]).T

    msm_df.columns=["tic_value","free_energy","protein_name","mdl_index"]

    if errorbars:
        bayes_df = bayes_one_dim_free_energy(prj,prt,tic_index, lin_spaced_tic=lin_spaced_tic)
        df = pd.concat([msm_df, bayes_df])
        return df
    else:
        return msm_df

def plot_2d(prj, prt, tic_list, x_array=None, y_array=None, n_bins=100):
    #basic sanity tests

    assert(len(tic_list)==2)
    if x_array is None:
        lin_spaced_tic = global_tic_boundaries([prt],tic_list,n_bins)
        x_array = lin_spaced_tic[tic_list[0]]
        y_array = lin_spaced_tic[tic_list[1]]
    else:
        n_bins = len(x_array) - 1

    H_overall = np.zeros((n_bins,n_bins))
   #get data
    c_x = prt.tic_dict[tic_list[0]]
    c_y = prt.tic_dict[tic_list[1]]
    for i in range(prt.n_states_):
        H,x,y=np.histogram2d(c_x[i],c_y[i],bins=[x_array,y_array],normed=True)
        H_overall = H_overall + prt.msm.populations_[i]*H
    H_copy = -0.6*np.log(H_overall)

    return H_copy.T

