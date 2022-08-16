#!/bin/bash

# calculates the rmsd for all complexes supplied 
# min rmsd script should be in the same directory as this wrapper
# this script should be in the same directory as all the complexes that dock6 was run on

for complex in $(<$1)
  do
  echo ${complex}
  #cd ${complex}
  file1=ligand_pdbqts/${complex}_ligand.pdbqt
  echo ${file1}
  file2=vina_modified_pdbqts/${complex}_ligand1_vina_out_mod.pdbqt
  echo ${file2}   
  ./min_rmsd.ksh $file1 $file2 > pdb_vina_results_pdbqt_proper_mod/${complex}.txt &
  #cd ..
done 
