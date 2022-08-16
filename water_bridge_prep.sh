#!/bin/bash

# prepare proteins and ligands for water bridge calculations. 

for dir in $(<$1)
do
cd ${dir}
obabel -imol2 ${dir}_ligand.mol2 -opdb -O ${dir}_ligand.pdb
grep "^ATOM" ${dir}_protein.pdb > ${dir}_protein_ligand.pdb
cat ${dir}_ligand.pdb | sed 's/ATOM/LIGATOM/' > ${dir}_ligand_mod.pdb
grep "^LIGATOM" ${dir}_ligand_mod.pdb >> ${dir}_protein_ligand.pdb
echo "TER" >> ${dir}_protein_ligand.pdb
grep "^HETATM*" ${dir}_protein.pdb >> ${dir}_protein_ligand.pdb
echo "END" >> ${dir}_protein_ligand.pdb
cd ..
done
