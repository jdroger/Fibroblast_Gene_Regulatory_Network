"""
Functions for refining GRN as inferred from Arboreto
"""

import os
import pandas as pd
import numpy as np

# Library filtering functions
def getLibPath(start_directory, sub_directory, filter_extension=None):
    """
    Using top-level and subdirectory, returns 
    string of filename matching extension
    """
    for root, _, files in os.walk(start_directory+sub_directory):
        for file in files:
            if filter_extension is None or file.lower().endswith(filter_extension):
                return os.path.join(root, file)


def importLibraries(dirpath, libraryname, both=True):
    """
    Using filename of GRN, imports libraries of TF-target interactions.
    :dirpath:   name of top-level directory containing library folders.
    :filename:  name of library to use 
                    Library options: CHEA, TRANSFAC, ENCODE
    """
    libextension = "_gene_attribute_edges.txt.gz"
    libfile = getLibPath(dirpath, libraryname, filter_extension=libextension)
    library = pd.read_csv(libfile, 
                        index_col=None, header=0, 
                        low_memory=False, sep=r"\s+")
    if both:
        libname = "TRANSFACpredicted"
        libfile = getLibPath(dirpath, libname, filter_extension=libextension)
        library2 = pd.read_csv(libfile, 
                        index_col=None, header=0, 
                        low_memory=False, sep=r"\s+")
        library = pd.concat([library,library2], axis=0).drop_duplicates()
    return library.drop(index=0)

def filterWithLibrary(grn, library):
    """given a grn df and library df, finds edges matching
    the library and returns filtered grn with matching edges"""
    grn_pairs = grn["TF"]+"-"+grn["target"]
    lib_pairs = library["target"]+"-"+library["source"]
    inLibrary = grn_pairs.isin(lib_pairs)
    grn_filtered = grn.loc[inLibrary, :]
    return grn_filtered


# =========================
# Input-output filtering functions

def findOutputs(grn, *regs, indivregs=None):
    """Given regex strings for output gene families, 
    finds outputs included in full GRN"""
    if indivregs is not None:
        reg_list = "|".join((regs) + indivregs)
    else:
        reg_list = "|".join((regs))

    outputs = grn.loc[grn["target"].str.contains(reg_list, regex=True).values,"target"].unique()
    return outputs


def findInputsOrOutputs(grn,values,column):
    """
    Given a list of gene names, 
    finds names located in either TF or target columns in GRN
    """
    vals_any = []
    for val in values:
        val_any = grn[column].isin([val]).any()
        vals_any.append(val_any)
    vals_dict = dict(zip(values,vals_any))
    # filter dictionary for true items
    keys = []
    for item in vals_dict.items():
        if item[1] == True:
            keys.append(item[0])
    return keys


def filterInputsOrOutputs(grn,keys,column):
    """Given a list of gene names and the corresponding
    column, returns a filtered GRN containing only
    elements in the list."""
    for key in keys:
        grn_key = grn.loc[grn[column]==key,:]
        # remove 'index' column if necessary 
        # (in order to avoid issues with concatenation)
        if grn.columns.isin(["index"]).any():
            grn_key = grn_key.drop(columns=["index"])
        if key == keys[0]:
            grn_filt = grn_key
        else:
            grn_filt = pd.concat([grn_filt, grn_key], axis=0)
    return grn_filt


def filterInOutNetwork(grn_in, grn_out, grn_targets, include_intermediates=False):
    """Given GRNs filtered for inputs TFs only, outputs targets only,
    and those included in libraries, function returns a subnetwork
    containing edges leading from inputs to outputs via intermediate
    TFs, with optional inclusion of indermediate (TF-TF) edges"""
    tfs = grn_targets["TF"].unique()
    in_tf = grn_in.loc[grn_in["target"].isin(tfs), :]
    tf_out = grn_out        # assumed b/c of TF column
    if include_intermediates:
        # find tf-tf edges connected to "input-tf" or "tf-output" edges
        grn_tfs = grn_targets.loc[grn_targets["target"].isin(grn_targets["TF"].unique()), :]
        inclInput = grn_tfs["TF"].isin(in_tf["target"])
        inclOutput = grn_tfs["target"].isin(tf_out["TF"])
        tftf_in = grn_tfs.loc[inclInput, :]
        tftf_out = grn_tfs.loc[inclOutput, :]
        tftf = tftf_in.loc[tftf_in.index.isin(tftf_out.index), :]
        grn_filtered = pd.concat([in_tf, tftf, tf_out], axis=0)
    else:
        grn_filtered = pd.concat([in_tf, tf_out], axis=0)
    return grn_filtered


# =========================
# Simultaneous top-down + bottom-up DFS scripts

def testPathInv(path,grn,rules):
    """Test that a found path meets bottom-up search requirements"""
    meetsRules = []
    for item in range(len(path)):
        if item < len(path)-1:
            target = path[item]
            TF = path[item+1]
            toTest = grn.loc[(grn["target"]==target) & (grn["TF"]==TF),"importance"]
            neighbors_all = grn.loc[grn["target"]==target,:]
            quant = neighbors_all["importance"].quantile(q=rules[1])
            if ((toTest >= rules[0]).bool()) | ((toTest >= quant).bool()):
                meetsRules.extend([True])
            else:
                meetsRules.extend([False])
    return all(meetsRules)


def findPathsBoth(grn,start,output_keys,rules=[1,0.5]):
    """
    Given a network of pairwise TF-target interactions, a 
    starting TF, and a set of output genes, 1) uses a top-down 
    search algorithm to find paths between the input and outputs, 
    and 2) checks that found paths meets the same rules for a 
    bottom-up search.
    """
    stack = [(start,[start])]
    while stack:
        (vertex,path) = stack.pop()
        neighbors_all = grn.loc[grn["TF"]==vertex,:]
        quant = neighbors_all["importance"].quantile(q=rules[1])
        toKeep = ((neighbors_all["importance"]>rules[0]) | (neighbors_all["importance"]>quant))
        neighbors = neighbors_all.loc[toKeep,"target"]
        for neigh in neighbors:
            if neigh not in path:
                if neigh in output_keys:
                    meetsInvRules = testPathInv(list(reversed(path + [neigh])),grn,rules)
                    if meetsInvRules:
                        yield path + [neigh]
                else:
                    stack.append((neigh, path + [neigh]))


def findPathImps(grn,pathlist):
    """
    
    """
    tot_imp = []
    vars_imp = []
    for paths in pathlist:
        path_imp = []
        for item in range(len(paths)):
            if item < len(paths)-1:
                TF = paths[item]
                target = paths[item+1]
                criteria = ((grn["TF"] == TF) & (grn["target"] == target))
                imp = grn.loc[criteria, "importance"]
                if item == 0:
                    path_imp = [float(imp)]
                else:
                    path_imp = path_imp + [float(imp)]
        tot_imp.extend([np.sum(path_imp)])
        vars_imp.extend([np.std(path_imp)])
    paths_imp = pd.DataFrame(data=[pathlist, tot_imp, vars_imp], 
                             index=["path", "importance_total", "importance_sd"]).transpose()
    
    # add additional metrics + metadata
    paths_imp["importance_mean"] = paths_imp["importance_total"] / paths_imp["path"].str.len()
    paths_imp["importance_cv"] = paths_imp["importance_sd"] / paths_imp["importance_mean"]
    paths_str = []
    for _, series in paths_imp.iterrows():
        paths_str.extend([''.join('->'+str(i) for i in series["path"])[2:]])
    paths_imp["path_string"] = paths_str
    paths_imp["input"] = paths_imp["path"].str[0]
    paths_imp["output"] = paths_imp["path"].str[-1]
    paths_imp["TF"] = paths_imp["path"].str[-2]
    return paths_imp


def findPathRows(grn,pathlist):
    """Extract interaction pairs from paths"""
    idxs = []
    for paths in pathlist:
        for item in range(len(paths)):
            if item < len(paths)-1:
                TF = paths[item]
                target = paths[item+1]
                criteria = ((grn["TF"]==TF) & (grn["target"]==target))
                idx = grn.loc[criteria,:].index
                idxs.extend(idx)
    paths_idx = grn.loc[idxs,:].sort_values("importance", ascending=False).drop_duplicates()
    return paths_idx







# =========================
# Runtime function

def refineGRN(grn, 
            libraryname, 
            dir_path="D:\\Research\\Aim3\\data_TFdatabases\\", 
            lib_both=True, 
            output_regex=False):
    """
    Top-level function for GRN refinement.
    :grn:           n x 3 pandas dataframe containing "TF", "target" and "importance"
    :libraryname:   string containing TF-target database used for pruning
                        Current options:    "CHEA", "TRANSFACpredicted", 
                                            "TRANSFACcurated", "ENCODE"
    :dir_path:      string containing top-level directory for all databases,
                    libraries should be located in folders matching libraryname
    :lib_both:      Boolean determining 'both' argument in importLibraries fcn
    :output_regex:  Boolean determining whether regex should be used in 
                    choosing gene outputs (i.e. gene families)
    :grn_final:     n x 3 pandas datafram containing refined edges
    """
    # import df
    # db_path = "CHEA\\"
    # file_name = "CHEA_both_GRN_09012020_EdgeR.csv"
    # grn = pd.read_csv(dir_path+db_path+file_name, header=0, index_col=0)
    # grn = grn.reset_index()
    # print(grn.shape)

    # filter for edges contained in librar(ies)
    library = importLibraries(dir_path, libraryname, both=lib_both)
    grn_targets = filterWithLibrary(grn, library)
    print(grn_targets.shape)

    # filter for edges connected to desired inputs/outputs
    if output_regex:
        reg_col = r"^COL\d{1,2}A\d$"
        reg_mmp = r"^MMP\d"
        reg_timp = r"TIMP\d"
        reg_cts = r"CTS+[A-Z]"
        reg_tgfb = r"TGFB\d$"
        reg_thbs = r"THBS\d"
        reg_lox = r"^LOX"
        reg_others = ("SPP1","POSTN","^FN1","SPARC$",
                    "TNC","CTGF","SERPINE1","ACTA2")

        outputs = findOutputs(grn, 
                            reg_col, reg_mmp, reg_timp, reg_cts, 
                            reg_tgfb, reg_thbs, reg_lox, 
                            indivregs=reg_others)
    else:
        outputs = ["CTGF","FN1","ACTA2","TIMP1","TIMP2","SERPINE1","MMP12",
                "MMP14","MMP1","MMP2","MMP3","MMP8","MMP9","POSTN","COL1A1",
                "COL1A2","COL3A1","TNC","THBS4","SPP1"]
    
    inputs = ["STAT1","STAT3","JUN","FOS","NFKB1","RELA","CREB1","CREBBP",
            "SMAD3","MYC","NFATC1","NFATC3","SRF","TEAD2","TEAD4","YAP1","WWTR1"]
    
    input_keys = findInputsOrOutputs(grn_targets, inputs, "TF")
    output_keys = findInputsOrOutputs(grn_targets, outputs, "target")
    grn_outputs = filterInputsOrOutputs(grn_targets, output_keys, "target")
    grn_inputs = filterInputsOrOutputs(grn_targets, input_keys, "TF")
    paths_all_inout = filterInOutNetwork(grn_inputs, grn_outputs, grn_targets)
    print(paths_all_inout.shape)

    # find input-output paths via modified DFS algorithm
    paths_found_bothsearch = []
    for inp in input_keys:
        paths_found_bothsearch.extend(list(findPathsBoth(paths_all_inout, 
                                                        inp, output_keys, 
                                                        rules=[1,0.75])))
    paths_found_bothsearch_imp = findPathImps(paths_all_inout, 
                                            paths_found_bothsearch)
    grn_final = findPathRows(paths_all_inout.reset_index(), 
                            paths_found_bothsearch_imp["path"])
    print(grn_final.shape)

    return grn_final