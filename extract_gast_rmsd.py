#!/usr/local/bin/python3

# extracts the gasteiger RMSDs 

import pandas as pd
import os

df = pd.DataFrame(columns = ["PDB_ID", "Gast_RMSD"])

os.chdir("pdb_gast_results_pdbqt_proper")
for file in os.listdir():
    pdb = file[:-4]
    print(pdb)
    line_count = 0
    with open(file, "r") as infile:
        for line in infile:
            line_count += 1
            if line_count == 1:
                rmsd_line = line
                break 
    rmsd = rmsd_line[10:14]
    df.loc[len(df)] = [pdb, rmsd]

df = df.sort_values('PDB_ID', ascending=True)
print(df)
os.chdir("..")
df.to_csv("Gast_pdb_rmsd_pdbqt_proper.tsv",sep="\t",header=True)
