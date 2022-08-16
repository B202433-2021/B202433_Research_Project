#!/usr/local/bin/python3

# count the number of water bridges between protein and ligand

import MDAnalysis
import MDAnalysis.analysis.hydrogenbonds
import os 
import sys
import subprocess

def get_complex_list():

    complex_list_file = sys.argv[1]
    with open(complex_list_file, "r") as infile:
        complexes = infile.read()
        complex_list = complexes.split("\n")
        complex_list = complex_list[:-1]
        return complex_list

def perform_analysis(filename):
    u = MDAnalysis.Universe(filename)
    #lig_res_command = f'grep "LIGATOM" {filename} | head -1 | cut -c 21-24'
    #lig_res_full = subprocess.check_output(lig_res_command, shell=True).decode("utf-8")
    #lig_res = lig_res_full[0:3]
    #print(lig_res)    
    w = MDAnalysis.analysis.hydrogenbonds.WaterBridgeAnalysis(u, 'record_type ATOM', 'record_type LIGATOM', 'resname HOH')
    w.run()
    print(w.results.timeseries)

comp_list = get_complex_list()

for comp in comp_list: 
    os.chdir(comp)
    print(comp)
    file = str(comp) + "_protein_ligand.pdb"
    perform_analysis(file)  
    os.chdir("..")
