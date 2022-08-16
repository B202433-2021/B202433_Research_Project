#!/usr/bin/ksh

# clean directories then run vina 

for dir in $(<$1); do
print $dir; sleep 1
    cp autoAD_test_new.ksh ${dir}/
    cd ${dir} 
    rm -f rm -f *.log *.txt *tmp* *.pdbqt *.dlg *pocket_box* *ligand1* vina.dpf *protein_H* *trim* INSPH *.in sphgen_cluster.pdb t.pdb *_ligand.pdb *_protein_site.pdb
    # remove waters and cofactors from protein pdb file 
    grep "^ATOM " ${dir}_protein.pdb > ${dir}_protein_cleaned.pdb
    mv ${dir}_protein_cleaned.pdb ${dir}_protein.pdb
    print "proteins have had waters and cofactors removed. Running ./autoAD_test_new.ksh"
    ./autoAD_test_new.ksh -cqjbv ${dir}_ligand.sdf ${dir}_protein.pdb ${dir}_pocket.pdb &
    cd ..
done
