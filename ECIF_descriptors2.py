#!/usr/local/bin/python3

import ecif
import pandas as pd
import os
import sys

# ecif.py and PDB_Atom_Keys.csv taken from https://github.com/DIFACQUIM/ECIF
# calculates the counts of all ECIFs and ligand based descriptors for the complexes

# get complex list
complex_list_file = sys.argv[1]
with open(complex_list_file, "r") as infile:
    complexes = infile.read()
    complex_list = complexes.split("\n")
    complex_list = complex_list[:-1]

# create dataframe with headers as the PDB ID, possible ECIFs and the ligand descriptors
PDBHeader = []
PDBHeader.append("PDB_ID")
ECIFHeaders = [header.replace(';','') for header in ecif.PossibleECIF]
print(len(ecif.PossibleECIF))
LigDesHeaders = ecif.LigandDescriptors
headers = PDBHeader + ECIFHeaders + LigDesHeaders
df = pd.DataFrame(columns = headers)
#print(df)

# loop over all files/directories in the directory
for dir in complex_list:
    os.chdir(dir)
    # get ECIF data and ligand descriptors and add this as a new row to the dataframe
    if str(dir) + "_protein_cleaned.pdb" in os.listdir() and str(dir) + "_ligand.mol" in os.listdir():
      rec = str(dir) + "_protein_cleaned.pdb"
      lig = str(dir) + "_ligand.mol"
      pdb = []
      pdb.append(str(dir))
      print(os.getcwd(), rec, lig)
      ECIF_data = ecif.GetECIF(rec, lig, distance_cutoff=6.0)
      Descriptors = ecif.GetRDKitDescriptors(lig)
      all_data = pdb + ECIF_data + list(Descriptors)
      df.loc[len(df)] = all_data
      os.chdir("..")
      # print(df)
    else:
      print(f"Ligand mol and clean protein pdb files not found in directory {dir}")
      os.chdir("..") 
  
print(df)
df.to_csv("ECIF_LigDes.tsv",sep="\t",header=True)
