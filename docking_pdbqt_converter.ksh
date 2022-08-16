#!/usr/bin/ksh
source /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/bin/mglenv.sh > /dev/null 2>&1

# converts files to pdbqt format

#cd gast_mol2_files
#for file in *.mol2
#cd vina
#for file in *.pdbqt
for dir in $(<$1)
do
# cd ${dir}
#dir=${file:0:4}
#/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l ${dir}_ligand.mol2 -o ../ligand_pdbqts/${dir}_ligand.pdbqt &
#/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l ${dir}_ligand1_out_1_MOPAC_DOCKed_scored.mol2 -o ${dir}_ligand1_out_1_MOPAC_DOCKed_scored.pdbqt
#/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l ${dir}_ligand1_out_1_gast_DOCKed_scored.mol2 -o ../gast_pdbqts/${dir}_ligand1_out_1_gast_DOCKed_scored.pdbqt 
/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l vina/${dir}_ligand1_vina_out.pdbqt -o vina_modified_pdbqts/${dir}_ligand1_vina_out_mod.pdbqt 
echo "${dir}"
#cd ..
done
