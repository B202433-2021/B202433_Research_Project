#!/bin/bash

# counts the number of hetatms 

for dir in $(<$1)
do cd ${dir}
echo "${dir}"
pdb=${dir}_protein.pdb
grep "^HETATM" ${pdb} | awk '{FS=""; if(!($18 == "H" && $19 == "O" && $20 == "H")){print $0;}}'  
cd ..
done
