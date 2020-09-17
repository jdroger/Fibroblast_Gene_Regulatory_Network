"""
Scripts for inferring GRNs, including:
- Single inference
- k-fold cross validation
- Conversion to Netflux models
"""

import os
# import time
import pandas as pd
import numpy as np
import asyncio

from distributed import Client, LocalCluster
from arboreto.algo import grnboost2

from src.GRNrefinement import refineGRN
import src.GRNvalidation as gv

asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())




def importData(filepath):
    data = pd.read_csv(filepath, index_col=0, header=1).dropna()
    return data


def importTFs(dirpath, libraryname, both=True):
    """
    Imports list of TFs from specified top-level directory
    (dirpath) and library string (libraryname). Optionally
    adds a second library to create composite list for 
    wider coverage of TFs.
    """
    from src.GRNrefinement import getLibPath
    libextension = "_attribute_list_entries.txt.gz"
    libfile = getLibPath(dirpath, libraryname, filter_extension=libextension)
    tf_all = pd.read_table(libfile)
    if both:
        libname = "TRANSFACpredicted"
        libfile = getLibPath(dirpath, libname, filter_extension=libextension)
        tf_2 = pd.read_table(libfile)
        tf_all = pd.concat([tf_all, tf_2], axis=0).drop_duplicates()
    return tf_all


def processData(data, cutoff=1, num_samples=2):
    """Given a pandas dataframe, returns a transposed numpy array
    of values and associated gene names for input into grnboost2. 
    Additionally applies a threshold to data as used for EdgeR.
    
    Note: assumes data is CPM data for thresholding purposes"""
    data_t = data.transpose()
    # apply threshold via EdgeR method
    data_threshold = data_t >= cutoff
    data_keep = data_t.loc[:, (data_threshold.sum(axis=0) >= num_samples).values]
    data_array = data_keep.values
    data_genes = data_keep.columns.values
    return data_array, data_genes


def saveGRN(grn, savedir, fold=None, 
            trainingset=None, testingset=None, 
            suffix=None):
    if suffix is not None:
        ext_csv = "_"+suffix+".csv"
        ext_txt = "_"+suffix+".txt"
    else:
        ext_csv = ".csv"
        ext_txt = ".txt"
    if fold is not None:
        filename_grn = "GRN_CV_fold"+str(fold)+ext_csv
        if trainingset and testingset is not None:
            filename_sets = "datasets_CV_fold"+str(fold)+ext_txt
            with open(savedir+filename_sets, "w") as output:
                output.write("_Training_\n")
                for train in trainingset[fold]:
                    output.write(train + "\n")
                output.write("_Testing_\n")
                for test in testingset[fold]:
                    output.write(test + "\n")
    else:
        filename_grn = "GRN_single"+ext_csv
    grn.to_csv(savedir+filename_grn)
    return



# runtime scripts
def inferGRN(filename, 
            libpath, libname, lib_both=True,
            savedir=None, suffix=None, seed=None):
    """
    Top-level script for inferring gene regulatory network
    from a given dataset using the Arboreto GRNboost2 algorithm.
    :filename:  path to CSV file containing gene expression data.
    :libpath:   path to directory containing sub-folders for TF-target
                libraries.
    :libname:   string of TF-target library used for inference.
    :lib_both:  (optional) Boolean operator determining use of additional
                library (TRANSFACpredicted) for wider TF coverage
    :savedir:   (optional) path to directory for saving final CSV.
    :seed:      (optional) integer for inference algorithm seed
    """

    # import cpm + library data
    cpm = importData(filename)
    cpm_array, cpm_genes = processData(cpm)

    tf_all = importTFs(libpath, libname, lib_both)
    tf_names = tf_all["GeneSym"].to_list()


    # setup Dask cluster
    client = Client(LocalCluster())
    print(client.dashboard_link)
    


    # infer + refine GRN
    grn = grnboost2(expression_data=cpm_array,
                    gene_names=cpm_genes,
                    tf_names=tf_names,
                    client_or_address=client,
                    seed=seed)
    grn_refined = refineGRN(grn, libname, dir_path=libpath)

    if savedir is not None:
        saveGRN(grn_refined, savedir, suffix=suffix)

    client.shutdown()
    return grn_refined


def crossvalidateGRN(filename, 
                    libpath, libname, k, lib_both=True,
                    savedir=None, suffix=None, seed=None):
    """
    Top-level script for k-fold cross validation of gene regulatory 
    network inference using the Arboreto GRNboost2 algorithm.
    :filename:  path to CSV file containing gene expression data.
    :libpath:   path to directory containing sub-folders for TF-target
                libraries.
    :libname:   string of TF-target library used for inference.
    :k:         integer specifying number of folds for CV
    :lib_both:  (optional) Boolean operator determining use of additional
                library (TRANSFACpredicted) for wider TF coverage
    :savedir:   (optional) path to directory for saving final CSV.
    :seed:      (optional) integer for inference algorithm seed
    """

    # import cpm + library data
    cpm = importData(filename)

    tf_all = importTFs(libpath, libname, lib_both)
    tf_names = tf_all["GeneSym"].to_list()

    # create and assign CV folds
    folds = gv.makeFolds(cpm, k)
    training, testing = gv.assignFolds(folds)

    # setup Dask cluster
    client = Client(LocalCluster())
    print(client.dashboard_link)

    # infer + refine GRN for each fold
    fold = 0
    while fold < k:
        cpm_fold = cpm.loc[:, training[fold]]
        cpm_array, cpm_genes = processData(cpm_fold)

        grn = grnboost2(expression_data=cpm_array,
                        gene_names=cpm_genes,
                        tf_names=tf_names,
                        client_or_address=client,
                        seed=seed)
        grn_refined = refineGRN(grn, libname, dir_path=libpath)

        if savedir is not None:
            saveGRN(grn_refined, savedir, fold=fold, suffix=suffix, 
                    trainingset=training, testingset=testing)
        
        # store all refined GRNs
        if fold == 0:
            grn_all = grn_refined
        else:
            grn_all = grn_all.append(grn_refined)

        fold = fold + 1

    client.shutdown()
    return grn_all