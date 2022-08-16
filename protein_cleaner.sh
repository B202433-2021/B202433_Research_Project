#!/bin/bash

# removes waters and all cofactors from a protein 

for dir in $(<$1)
do 
cd ${dir}
grep "^ATOM " ${dir}_protein.pdb > ${dir}_protein_cleaned.pdb &
cd ..
done
