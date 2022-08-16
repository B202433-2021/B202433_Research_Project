#!/usr/bin/ksh
source /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/bin/mglenv.sh > /dev/null 2>&1

# converts crystallographic and MOPAC docked ligand to pdbqt. Also converts protein to pdbqt.

function prep_lig {
  # converts ligand to pdbqt but keeps non polar hydrogens
  /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l $1_ligand.mol2 -U lps -o $1_ligand.pdbqt 
}

function prep_protein {
  # adds position optimised hydrogens to the protein then converts to pdbqt. Keeps non polar hydrogens
  grep -v "^TER" $1_protein_cleaned.pdb | grep -v "^REMARK" | grep "^[AH][TE][OT][MA][ T][ M]" | sed "s/^\(.\{12\}\)SE   MSE/\1 S   MET/" | sed "s/^\(.\{17\}\)MSE/\1MET/" | cut -c1-76 > $1_protein_tmp.$$
  
  OLDPYTHONHOME=$PYTHONHOME
  unset PYTHONHOME
  print "/usr/bin/python /localdisk/home/s2173344/Research_Project/pdb2pqr-src-2.1.1/pdb2pqr.py  --ff=AMBER --with-ph=7.4  $1_protein_tmp.$$ $1_protein_H.pqr"
  /usr/bin/python /localdisk/home/s2173344/Research_Project/pdb2pqr-src-2.1.1/pdb2pqr.py  --ff=AMBER --with-ph=7.4  $1_protein_tmp.$$ $1_protein_H.pqr
  export PYTHONHOME=$OLDPYTHONHOME

  if [[ ! -s $1_protein_H.pqr ]]; then
  print "WARNING: PDB2PQR has failed, switching to OpenBabel..."
  print $1 > ../pdb2pqr_failures.txt
  print obabel -ipdb $1_protein_tmp.$$ -opdb -O $1_protein_H.pqr -p 7.4
  obabel -ipdb $1_protein_tmp.$$ -opdb -O $1_protein_H.pqr -p 7.4
  fi
  rm -f $1_protein_tmp.$$
  cut -c1-54 $1_protein_H.pqr > $1_protein_H.pdb
  print "/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -U lps -r $1_protein_H.pdb"
  /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -U lps -r $1_protein_H.pdb
  grep -v '^ATOM.*  H   .*0.000 HD$'  ${dir}_protein_H.pdbqt  > tmp.pdbqt.$$; mv tmp.pdbqt.$$ ${dir}_protein_H.pdbqt
}

for dir in $(<$1)
do
cd ${dir}
echo ${dir}
prep_lig ${dir} &
prep_protein ${dir} & 
cd ..
done
