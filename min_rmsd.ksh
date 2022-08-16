#!/usr/bin/ksh
#Gives RMSD and min. RMSD for two pdb or pdbqt files. Atoms have to be in the same order.
#Only works for pdb or pdbqt files - rmsd.ksh file1.pdb file2.pdb

function errhan {
print "USAGE: rmsd_min.ksh file1.pdbqt file2.pdbqt
Only works for pdb or pdbqt files that contain Autodock atom typing.
The first figure output is the minimum RMSD taking into account symmetry; atoms can be in any order.
The second figure output is the 'strict' RMSD; for this value to be correct the atoms must be in the same order in both files."
print $error
exit 1
}

if [[ ${#} == 0 ]]; then
  error="No input files detected"
  errhan
fi

if [[ $(grep -s "REMARK VINA RESULT:" $1) ]]; then
  print "Vina file format detected for $1;"
  print "Assuming you want to compare the first docking pose listed in the file ..."
  block1=$(egrep "^[AH][TE][OT][MA][ T][ M]|^MODEL|^ENDMDL" $1  | awk '/'MODEL\ 1'$/,/ENDMDL/' | grep "^[AH][TE][OT][MA][ T][ M]")
else
  block1=$(grep  "^[AH][TE][OT][MA][ T][ M]" $1 )
fi
if [[ $(grep -s "REMARK VINA RESULT:" $2) ]]; then
  print "Vina file format detected for $2;"
  print "Assuming you want to compare the first docking pose listed in the file ..."
  block2=$(egrep "^[AH][TE][OT][MA][ T][ M]|^MODEL|^ENDMDL" $2  | awk '/'MODEL\ 1'$/,/ENDMDL/' | grep "^[AH][TE][OT][MA][ T][ M]")
else
  block2=$(grep  "^[AH][TE][OT][MA][ T][ M]" $2 )
fi


types1=($(print "$block1" | cut -c 78-80 ))
xcoords1=($(print "$block1" | cut -c31-38 ))
ycoords1=($(print "$block1" | cut -c39-46 ))
zcoords1=($(print "$block1" | cut -c47-54 ))
#print number of variables is ${#xcoords1[*]} and last variable is ${xcoords1[$((${#xcoords1[*]}-1))]}

types2=($(print "$block2" | cut -c 78-80 ))
xcoords2=($(print "$block2" | cut -c31-38 ))
ycoords2=($(print "$block2" | cut -c39-46 ))
zcoords2=($(print "$block2" | cut -c47-54 ))
#print number of variables is ${#xcoords2[*]} and last variable is ${xcoords2[$((${#xcoords1[*]}-1))]}

#fix for pdbqts with NA instead of N
for ((c=0; c < ${#xcoords1[*]}; c++)); do
  if [[ ${types1[$c]// } == NA ]]; then
    types1[$c]=N
  fi
  
  if [[ ${types2[$c]// } == NA ]]; then
    types2[$c]=N
  fi
  
  if [[ ${types1[$c]// } == SA ]]; then
    types1[$c]=S
  fi
  
  if [[ ${types2[$c]// } == SA ]]; then
    types2[$c]=S
  fi
  
done

shortest_distance=100000000
#Scan all atoms in the reference structure and select the nearest atom of identical type to that in the query structure
for ((c=0; c < ${#xcoords1[*]}; c++)); do
  if [[ ${types1[$c]// } != HD ]]; then
    unset distance
    shortest_distance=1000000
    for ((n=0; n < ${#xcoords1[*]}; n++)); do
#         print ${types1[$c]} == ${types2[$n]}
      if [[ ${types1[$c]// } == ${types2[$n]// } ]]; then
        ((distance=((${xcoords1[$c]}-(${xcoords2[$n]}))**2)+((${ycoords1[$c]}-(${ycoords2[$n]}))**2)+((${zcoords1[$c]}-(${zcoords2[$n]}))**2) ))
#        print $distance $shortest_distance
        if (( distance < shortest_distance )); then
          shortest_distance=$distance
          type1=$c
          type2=$n
        fi
      fi
    done
  fi
#  print shortest distance between types ${types1[$type1]} and ${types2[$type2]} is $shortest_distance
  #print $shortest_distance
  ((min_sum+=$shortest_distance))
done
typeset -F2 min_rmsd
((min_rmsd=sqrt( (1/${#xcoords1[*]}.0)*$min_sum) ))


#for ((c=0; c < ${#xcoords1[*]}; c++)); do
#  print "variance=((${xcoords1[$c]}-(${xcoords2[$c]}))**2)+((${ycoords1[$c]}-(${ycoords2[$c]}))**2)+(($#{zcoords1[$c]}-(${zcoords2[$c]}))**2 )))"
#  ((variance=((${xcoords1[$c]}-(${xcoords2[$c]}))**2)+((${ycoords1[$c]}-(${ycoords2[$c]}))**2)+(($#{zcoords1[$c]}-(${zcoords2[$c]}))**2 )))
#  ((sum+=$variance))
#done
#print $sum ${#xcoords1[*]} 

#typeset -F2 rmsd
#((rmsd=sqrt( (1/${#xcoords1[*]}.0)*$sum) ))
print "Min_RMSD: $min_rmsd Strict_RMSD: rmsd"
