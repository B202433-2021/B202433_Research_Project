#!/bin/bash

# converts mol2 files to mol format 

for dir in $(<$1)
do 
cd ${dir}
obabel -imol2 ${dir}_ligand.mol2 -omol -O ${dir}_ligand.mol &
cd ..
done
