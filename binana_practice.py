#!/usr/local/bin/python3

import sys  
sys.path.append('/localdisk/home/s2173344/Research_Project/binana-2.1/python')
import binana
import pandas as pd
import os

headers = []

complex_list_file = sys.argv[1]
with open(complex_list_file, "r") as infile:
    complexes = infile.read()
    complex_list = complexes.split("\n")
    complex_list = complex_list[:-1]

# generate every possible interaction type present in the whole set of samples - put in the "headers" list
for complex in complex_list:
    os.chdir(complex)
    lig_file = complex + "_ligand.pdbqt"
    rec_file = complex + "_protein.pdbqt"
    ligand, receptor = binana.load_ligand_receptor.from_files(lig_file, rec_file)
    all_inf = binana.interactions.get_all_interactions(ligand, receptor)
    for int_type in all_inf.keys():
        # handle the error when there are no interactions of that type
        try:
            for int_subtype in all_inf[int_type]['counts']:
                header = int_type + "_" + int_subtype
                headers.append(header)
        except TypeError:
            continue
    os.chdir("..")    

#print("Headers are: ", headers, "\n")
uniq_headers = sorted(set(headers))
#print(uniq_headers)

uniq_headers.insert(0, "PDB_ID")
df = pd.DataFrame(columns = uniq_headers)
pdb_id = "1atw"
print(df)

for complex in complex_list:
    os.chdir(complex)
    count_row = []
    count_row.append(complex)
    
    lig_file = complex + "_ligand.pdbqt"
    rec_file = complex + "_protein.pdbqt"
    ligand, receptor = binana.load_ligand_receptor.from_files(lig_file, rec_file)
    all_inf = binana.interactions.get_all_interactions(ligand, receptor)

    # get all the interactions for that complex in a dictionary
    interaction_dict = {}
    for int_type in all_inf.keys():
        # handle the error when there are no interactions of that type
        try:
            for int_subtype in all_inf[int_type]['counts']:
                header = int_type + "_" + int_subtype
                count = all_inf[int_type]['counts'][int_subtype]
                interaction_dict[header] = count
        except TypeError:
            continue
    print("Interaction dictionary for that sample is: ", interaction_dict, "\n")

    # loop over all possible interactions. Add the count if that interaction is present for that complex otherwise put 0
    for interaction in uniq_headers[1:]:
        if interaction in interaction_dict.keys():
            count_row.append(int(interaction_dict[interaction]))
        else:
            count_row.append(0)
    print(len(headers))
    print(len(count_row))
    print("Count row is: ", count_row)
    df.loc[len(df)] = count_row
    os.chdir("..")

print(df)
df.to_csv("binana.tsv",sep="\t",header=True)

