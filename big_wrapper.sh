#!/bin/bash

# script to run mopac docking in chunks to avoid server overload

for i in {50..1800..50}
do
start=${i}
echo "Complexes from ${start} to ${i}"  
cat comps_for_mopac.txt | head -$i | tail -50 > complexes.txt
cat complexes.txt
./vina_wrapper.sh complexes.txt
# ./vina_pdb_convert_wrapper.ksh complexes.txt
#./dock6.ksh complexes.txt
done
