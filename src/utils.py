"""
Utility functions used in inference: saveNetflux
"""

import pandas as pd

def saveNetflux(grn, savefile, composite=False, suffix=None): 
    # format reactions df
    grn = grn.reset_index()
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