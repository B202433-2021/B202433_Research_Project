#!/usr/local/bin/python3

# extracts the mopac RMSD 

import pandas as pd
import os

df = pd.DataFrame(columns = ["PDB_ID", "MOPAC_RMSD"])

os.chdir("pdb_MOPAC_results")
for file in os.listdir():
    pdb = file[:-4]
    line_count = 0
    with open(file, "r") as infile:
        for line in infile:
            line_count += 1
            if line_count == 2:
                rmsd_line = line
                break 
    rmsd = rmsd_line[10:14]
    df.loc[len(df)] = [pdb, rmsd]

print(df)
os.chdir("..")
df.to_csv("MOPAC_pdb_rmsd.tsv",sep="\t",header=True)
