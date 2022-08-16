#!/usr/local/bin/python3

# identifies the ECIF atom types present in the set of complexes 

import ecif
import pandas as pd
import os
import sys
from rdkit import Chem
from scipy.spatial.distance import cdist
from itertools import product
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator

# ecif.py and PDB_Atom_Keys.csv taken from https://github.com/DIFACQUIM/ECIF

def get_complex_list():

    complex_list_file = sys.argv[1]
    with open(complex_list_file, "r") as infile:
        complexes = infile.read()
        complex_list = complexes.split("\n")
        complex_list = complex_list[:-1]
        return complex_list

def GetAtomType(atom):
# This function takes an atom in a molecule and returns its type as defined for ECIF

    AtomType = [atom.GetSymbol(),
                str(atom.GetExplicitValence()),
                str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() != "H"])),
                str(len([x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "H"])),
                str(int(atom.GetIsAromatic())),
                str(int(atom.IsInRing())),
               ]

    return(";".join(AtomType))


def LoadSDFasDF(SDF):
# This function takes an SDF for a ligand as input and returns it as a pandas DataFrame with its atom types labeled according to ECIF

    m = Chem.MolFromMolFile(SDF, sanitize=False)
    m.UpdatePropertyCache(strict=False)
    at_types = []

    for atom in m.GetAtoms():
        if atom.GetSymbol() != "H": # Include only non-hydrogen atoms
            at_types.append(GetAtomType(atom))
    return at_types
  
comp_list = get_complex_list()

all_atom_types = []

for comp in comp_list:
    os.chdir(comp)
    lig = str(comp) + "_ligand.mol"
    print(lig)
    atom_types = LoadSDFasDF(lig)
    all_atom_types += atom_types
#    print(all_atom_types)
    os.chdir("..")

possible_atom_types = sorted(set(all_atom_types))
print(possible_atom_types)
print(len(possible_atom_types))
