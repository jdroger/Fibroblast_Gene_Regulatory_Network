"""
Utility functions used in inference: saveNetflux
"""

import pandas as pd

def saveNetflux(grn, savefile, inputs=None, suffix=None): 
    # format reactions df
    grn = grn.reset_index()

    if inputs:
        inputrxns = pd.DataFrame(data=" => " + inputs, columns=["Rule"])
        inputrxns["module"] = "input"
        inputrxns["ID"] = "i"+str(inputrxns.index)
        inputrxns["Weight"] = 1
        inputrxns["n"] = 1.4
        inputrxns["EC50"] = 0.6
        inputrxns["source"] = "GRN-inputs"
        inputrxns["notes"] = ""

    grn["module"] = "txn"
    grn["ID"] = "r" + str(grn.index)
    grn["Rule"] = grn["TF"] + " => " + grn["target"]
    grn["Weight"] = 1
    grn["n"] = 1.4
    grn["EC50"] = 0.6
    grn["source"] = "GRN"
    grn["notes"] = ""
    reactions = grn[["module", "ID", "Rule", "Weight", 
                    "n", "EC50", "source", "notes"]]
    
    if inputs:
        reactions = pd.concat([inputrxns, reactions], axis=0)
    
    # create species df from unique TF/target nodes
    speciesnames = pd.concat([grn["TF"], grn["target"]], 
                        axis=0, ignore_index=True).drop_duplicates().sort_values()
    species = pd.DataFrame(speciesnames, columns=["ID"])
    species["module"] = "txn"
    species["name"] = ""
    species["Yinit"] = 0
    species["Ymax"] = 1
    species["tau"] = 0.1
    species["type"] = "mRNA"
    species = species[["module", "ID", "name", "Yinit", "Ymax", "tau", "type"]]
    
    # export to xlsx
    with pd.ExcelWriter(savefile) as writer:
        species.to_excel(writer, sheet_name='species', index=False, startrow=1)
        reactions.to_excel(writer, sheet_name='reactions', index=False, startrow=1)

    return