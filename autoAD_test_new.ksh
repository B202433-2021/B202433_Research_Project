#!/usr/bin/ksh
source /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/bin/mglenv.sh > /dev/null 2>&1 

# performs the Vina preparation. Called like ./autoAD_test_new.ksh -cqjbv in this study. Vina is run through autoAD.ksh

USAGE=$'[-?\n@(#)$Id: autoAD'
USAGE+=$' 1.2.1  $\n]'
USAGE+="[-author?Douglas R. Houston <dhouston@staffmail.ed.ac.uk>]"
USAGE+="[-copyright?Copyright (c) D. R. Houston 2009.]"
USAGE+="[+NAME?autoAD --- Autodock automator]" 
USAGE+="[+DESCRIPTION?Prepares an SD file of compounds and "
USAGE+="a protein receptor for a Vina or Autodock run, can run either "
USAGE+="if option is specified.]"
USAGE+="[+EXAMPLE?autoAD.ksh -cqpgmds sdffile.sdf receptor.pdb sitepoints.mol2]"
USAGE+="[c:convert?Uses Babel to convert specified sdf file to individual mol2s.]"
USAGE+="[q:pdbqt?Uses AutodockTools to convert specified mol2 file to pdbqt.]"
USAGE+="[a:autodock?Run Autodock using specified dpf file.]"
USAGE+="[l:listdock?Run Autodock using specified docking list file.]"
USAGE+="[y:dyndock?Run Autodock using dynamically generated list (Dynamic mode).]"
USAGE+="[d:makedpfs?Make Autodock docking parameter files for every compound.]"
USAGE+="[j:vinadpf?Make Vina docking parameter file.]"
USAGE+="[l:prepose?Do not randomise starting pose.]"
USAGE+="[p:preprecep?Hydrogenate receptor using pdb2pqr, add charges and merge nonpolar hydrogens.]"
USAGE+="[g:makegpf?Make Autodock grid parameter file for receptor.]"
USAGE+="[m:makemaps?Calculate Autodock grids.]"
USAGE+="[x:maxparams?Set parameters to maximum (default is HTS settings).]"
USAGE+="[w:modes?Increase the maximum number of Vina binding modes output to 20 (default is 9).]"
USAGE+="[s:autodist?Job farm autodist in dynamic mode with Autodock.]"
USAGE+="[v:vinadist?Job farm Vina using dynamically generated list (Dynamic mode).]"
USAGE+="[h:listvina?Run Autodock using specified docking list file.]"
USAGE+="[z:dynvina?Run autodist in dynamic mode with Vina.]"
USAGE+="[n:vina?Run vina using specified ligand, receptor and binding site definition.]"
USAGE+="[e:evryting?Activate all preparatory options, run Autodock in Dynamic mode.]"
USAGE+="[b:vevryting?Activate all preparatory options, run Vina in Dynamic mode.]"
USAGE+="[o:prepDOCK6?Convert protein and compoundlibrary.sdf to individual .mol2 files for DOCK, calculate energy grids.]"
USAGE+="[u:makeDOCK6dpf?Make DOCK docking parameter files for every compound.]"
USAGE+="[f:dynDOCK6?Run autodist in dynamic mode with DOCK 6.8.]"
USAGE+="[k:distDOCK6?Job farm DOCK 6.9 using dynamically generated list (Dynamic mode). Takes no argument; assumes -ou has ran. Does not work on itioc1.]"
USAGE+="[t:verbose?Activate verbose mode (for bug fixing).]"
USAGE+=$'\n\n <compoundlist.sdf> <receptor.pdb> <gridcentre.pdb>\n\n'

scriptname=$0
function errhan {
print
eval $scriptname --man
print "\n$error"
exit 1
}


if (( $# < 1 )) || [[ ${1:0:1} != "-" ]]; then 
  error="You must specify some arguments and option flags."
  errhan $0  
fi

while getopts "$USAGE" optchar ; do
  case $optchar in
     h) listvina=listvina ;;
     l) pose=1 ;;
     c) sdf2pdbs=sdf2pdbs ;;
     q) pdbs2pdbqts=pdbs2pdbqts ;;
     d) makedpfs=makedpfs ;;
     u) makeDOCK6dpf=makeDOCK6dpf ;;
     p) preprecep=preprecep ;;
     g) makegpf=makegpf ;;
     m) makemaps=makemaps ;;
     x) max_params="--exhaustiveness=32" ; print Maximum search parameters activated. ;;
     w) num_modes="--num_modes 20" ;;
     a) autodock=autodock ;;
     l) autodock=listdock ;;
     y) autodock=dyndock ;;
     s) autodock=distdock ;;
     v) vina=vina ;;
     z) dynvina=dynvina ;;
     n) vina_oneoff=vina_oneoff ;;
     j) vinadpf=vinadpf ;;
     e) sdf2pdbs=sdf2pdbs 
        pdbs2pdbqts=pdbs2pdbqts
        makedpfs=makedpfs
        preprecep=preprecep 
        makegpf=makegpf
        makemaps=makemaps 
        autodock=distdock ;;
     b) if (( $# < 3 )); then
          error="Not enough arguments"
          errhan
        fi
        sdf2pdbs=sdf2pdbs 
        pdbs2pdbqts=pdbs2pdbqts
        preprecep=preprecep 
        vinadpf=vinadpf
        vina=vina ;;
     o) if (( $# < 4 )); then
          error="Not enough arguments; need compoundlist.sdf, receptor.pdb and bindsitedef.mol2"
          errhan
        fi
        prepDOCK6=prepDOCK6 ;;
     f) if (( $# < 2 )); then
          error="Not enough arguments. Do not use manually (you want -k)."
          errhan
        fi
        autodock=dynDOCK6 ;;
     k) if (( $# < 1 )); then
          error="Not enough arguments. Needs library.sdf"
          errhan
        fi
        autodock=distDOCK6 ;;
     t) verbose=1 
        print "Verbose mode activated" ;;
     *) error="Unable to recognize option ${1:-your arguments}." 
        errhan $0 ;;
  esac
done
shift $(($OPTIND - 1))

for arg in $*; do
  if [[ ! -e $arg || ! -f $arg ]]; then
    error="Cannot find file $arg"
    errhan
  fi
  ((c++))
  arglist+=($arg)
done

spacing=0.2

ligid=${arglist[0]%.*}
sortkey=$(print ${ligid} | wc -c )


#Build hostnames file
cores="bioinfmsc7.bio.ed.ac.uk 128"

if [[ ! -f cores.txt ]]; then
  print "No cores.txt found; using: \n$cores"
  print "$cores" > cores.txt
fi
while read -r line; do
#  print .${line##* }.
  for ((c = 0; c < ${line##* }; c++)); do
    hostnames+="${line%% *} "
  done
done < cores.txt
print "$hostnames" > hostnames.txt



if [[ $(uname -s) != Linux ]]; then
  error="Only works on Linux, sorry."
  errhan $0 
fi


function sdf2pdbs {
#Use babel to convert sdf into multiple pdbs, strips H out first (to remove mispositioned H) then re-adds it
if [[ ${arglist[0]##*.} != sdf && ${arglist[0]##*.} != pdb ]] && [[ ${arglist[0]##*.} != mol2 ]]; then
  error="Please specify a list of compounds in SD file format. "
  errhan
fi
print Checking ${arglist[0]##*.} file for errors ...
dot=$(head -1 ${arglist[0]} | cut -f1 -d" ")
sed -i "s/^$dot/${dot#.}/" ${arglist[0]} 
sed -i "s/0999 V2000/0   1 V2000/" ${arglist[0]} 
if [[ -z $(head -10000 ${arglist[0]} | grep "M  END") ]]; then
  print "Adding \"M  END\" to SD file ..."
  sed -i '/>  <PARTIALQ_INFO>/i\
M  END' ${arglist[0]}
fi
print ${arglist[0]%.sdf} > ligid.txt
print Removing any hydrogen atoms from ${arglist[0]##*.} file ...
obabel -i${arglist[0]##*.} ${arglist[0]} -o${arglist[0]##*.} -O ${ligid}_noH_tmp.${arglist[0]##*.}.$$ -d
print Adding hydrogen atoms to compounds, and converting to separate mol2 files ...
obabel -i${arglist[0]##*.} ${ligid}_noH_tmp.${arglist[0]##*.}.$$ -omol2 -O ${ligid}.mol2 -p 7.4 -m
#print Converting ${ligid}_OBH.sdf to individual mol2 files ...
#/usr/people/douglas/programs/openbabel/bin/babel -isdf ${ligid}_OBH.sdf -omol2 ${ligid}.mol2 -m
#print Adding hydrogens to mol2 files
#/usr/people/douglas/programs/openbabel/bin/babel -h -imol2 tmp$$*.mol2 -omol2 -m ${ligid}.mol2
#for file in $(ls ${ligid}tmp$$*.mol2); do
#  mv $file ${file/tmp$$}
#done
#rm *tmp*$$*
}


function pdbs2pdbqts {
#Use Autodock script to convert mol2s into pdbqts
if [[ -n $(ls -1 ${ligid}+([[:digit:]]*).mol2 2> /dev/null) ]]; then
  for file in $( ls -1 ${ligid}+([[:digit:]]*).mol2 | sort -n --key=1.$sortkey); do
    print "Converting $file to pdbqt format ..."
    /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l $file  > /dev/null 2>&1; wait
    if [[ ! -e ${file%.mol2}.pdbqt ]]; then
      print prepare_ligand4.py failed on $file, attempting to use PDB format instead ...
      obabel -imol2 $file -opdb -O ${file%.mol2}.pdb > /dev/null 2>&1 
      /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l ${file%.mol2}.pdb > /dev/null
    fi
  done
else
  for file in $( ls -1 ${ligid}?([[:digit:]]*).mol2 | sort -n --key=1.$sortkey); do
    print "Converting $file to pdbqt format ..."
    print "/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l $file"
    /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l $file  > /dev/null 2>&1; wait
    if [[ ! -e ${file%.mol2}.pdbqt ]]; then
      print prepare_ligand4.py failed on $file, attempting to use PDB format instead ...
      obabel -imol2 $file -opdb -O ${file%.mol2}.pdb > /dev/null 2>&1 
      /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -F -l ${file%.mol2}.pdb > /dev/null
    fi
  done
fi
print "Atom types in ligand files:"
cut -c78-79 ${ligid}*.pdbqt | sort -u
#for atomtype in $(cut -c78-79 ${ligid}*.pdbqt | sort -u); do
#  print "map ${arglist[1]%.*}_H.${atomtype}.map #atom-specific affinity map"
#done 
#print
}


function preprecep {
if [[ ${#arglist[*]} > 1 ]] && [[ ${arglist[1]##*.} == pdb ]]; then
  argno=1
elif [[ ${arglist[0]##*.} == pdb ]]; then
  argno=0
else
  error="Are you sure your arguments are in the right order?"
  errhan
fi
print "Preparing receptor ... "
grep -v "^TER" ${arglist[${argno}]} | grep -v "^REMARK" | grep "^[AH][TE][OT][MA][ T][ M]" | sed "s/^\(.\{12\}\)SE   MSE/\1 S   MET/" | sed "s/^\(.\{17\}\)MSE/\1MET/" | cut -c1-76 > ${arglist[${argno}]%.*}_tmp.$$

OLDPYTHONHOME=$PYTHONHOME
unset PYTHONHOME
print "/usr/bin/python /localdisk/home/s2173344/Research_Project/pdb2pqr-src-2.1.1/pdb2pqr.py  --ff=AMBER --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr"
/usr/bin/python /localdisk/home/s2173344/Research_Project/pdb2pqr-src-2.1.1/pdb2pqr.py  --ff=AMBER --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr

#export PYTHONHOME=/usr/people/douglas/programs/Python-2.7.11/build/
#/usr/people/douglas/programs/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py --ff=PARSE --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr
export PYTHONHOME=$OLDPYTHONHOME

if [[ ! -s ${arglist[${argno}]%.*}_H.pqr ]]; then
  print "WARNING: PDB2PQR has failed, switching to OpenBabel..."
  print obabel -ipdb ${arglist[${argno}]%.*}_tmp.$$ -opdb -O ${arglist[${argno}]%.*}_H.pqr -p 7.4
  obabel -ipdb ${arglist[${argno}]%.*}_tmp.$$ -opdb -O ${arglist[${argno}]%.*}_H.pqr -p 7.4
fi
rm -f ${arglist[${argno}]%.*}_tmp.$$
cut -c1-54 ${arglist[${argno}]%.*}_H.pqr > ${arglist[${argno}]%.*}_H.pdb
print "/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -U nphs_lps -r ${arglist[${argno}]%.*}_H.pdb"
/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -U nphs_lps -r ${arglist[${argno}]%.*}_H.pdb
grep -v '^ATOM.*  H   .*0.000 HD$'  ${arglist[${argno}]%.*}_H.pdbqt  > tmp.pdbqt.$$; mv tmp.pdbqt.$$ ${arglist[${argno}]%.*}_H.pdbqt
}


function boxsize {
#USAGE boxsize ligand spacing
if [[ ${1##*.} == mol2 ]]; then
  xcoords=($(awk '/^@<TRIPOS>ATOM/,/@<TRIPOS>BOND/' $1 | cut -c17-26 | sort -g))
  ycoords=($(awk '/^@<TRIPOS>ATOM/,/@<TRIPOS>BOND/' $1 | cut -c27-36 | sort -g))
  zcoords=($(awk '/^@<TRIPOS>ATOM/,/@<TRIPOS>BOND/' $1 | cut -c37-46 | sort -g))
# print number of variables is ${#xcoords[*]} and last variable is ${xcoords[$((${#xcoords[*]}-1))]} 1>&2
  pad=8
else
  xcoords=($(grep  "^[AH][TE][OT][MA][ T][ M]" $1 | cut -c31-38 | sort -g))
  ycoords=($(grep  "^[AH][TE][OT][MA][ T][ M]" $1 | cut -c39-46 | sort -g))
  zcoords=($(grep  "^[AH][TE][OT][MA][ T][ M]" $1 | cut -c47-54 | sort -g))
# print number of variables is ${#ycoords[*]}, first and last variables are ${ycoords[0]} and ${ycoords[$((${#ycoords[*]}-1))]} 1>&2
  pad=12
fi

((centrex=(${xcoords[0]}+${xcoords[$((${#xcoords[*]}-1))]})/2))
((maxx=${xcoords[$((${#xcoords[*]}-1))]}-(${xcoords[0]})))
#print maxx is $maxx 1>&2

((centrey=(${ycoords[0]}+${ycoords[$((${#ycoords[*]}-1))]})/2))
((maxy=${ycoords[$((${#ycoords[*]}-1))]}-(${ycoords[0]})))
#print maxy=${ycoords[$((${#ycoords[*]}-1))]}-(${ycoords[0]})  1>&2
#print maxy is $maxy 1>&2

((centrez=(${zcoords[0]}+${zcoords[$((${#zcoords[*]}-1))]})/2))
((maxz=${zcoords[$((${#zcoords[*]}-1))]}-(${zcoords[0]})))

#print Centre of the molecule is $centrex $centrey $centrez

typeset -i halfnpntsx
((halfnpntsx=((maxx+pad)/$2)/2))
((npntsx=halfnpntsx*2))

typeset -i halfnpntsy
((halfnpntsy=((maxy+pad)/$2)/2))
((npntsy=halfnpntsy*2))

typeset -i halfnpntsz
((halfnpntsz=((maxz+pad)/$2)/2))
((npntsz=halfnpntsz*2))

((space=($npntsx*$2)*($npntsy*$2)*($npntsz*$2) ))
if (( space > 47000 )); then
  print "WARNING: The search space volume is $space; > 47000 Angstrom^3!" 1>&2
  sleep 5
else
  print "All is well; the search space volume is $space; <= 47000 Angstrom^3" 1>&2
fi

print "npts $npntsx $npntsy $npntsz"
printf "gridcenter %3.3f %3.3f %3.3f\n" $centrex $centrey $centrez

print "HEADER    CORNERS OF BOX
REMARK    CENTER (X Y Z)  $centrex $centrey $centrez
REMARK    DIMENSIONS (X Y Z)   $(($npntsx*$2)) $(($npntsy*$2)) $(($npntsz*$2)) " > ${1%.*}_box.pdb

corners+=( $(($centrex-($halfnpntsx*$2))) $(($centrey-($halfnpntsy*$2))) $(($centrez-($halfnpntsz*$2))) )
corners+=( $(($centrex+($halfnpntsx*$2))) $(($centrey-($halfnpntsy*$2))) $(($centrez-($halfnpntsz*$2))) )
corners+=( $(($centrex+($halfnpntsx*$2))) $(($centrey-($halfnpntsy*$2))) $(($centrez+($halfnpntsz*$2))) )
corners+=( $(($centrex-($halfnpntsx*$2))) $(($centrey-($halfnpntsy*$2))) $(($centrez+($halfnpntsz*$2))) )
corners+=( $(($centrex-($halfnpntsx*$2))) $(($centrey+($halfnpntsy*$2))) $(($centrez-($halfnpntsz*$2))) )
corners+=( $(($centrex+($halfnpntsx*$2))) $(($centrey+($halfnpntsy*$2))) $(($centrez-($halfnpntsz*$2))) )
corners+=( $(($centrex+($halfnpntsx*$2))) $(($centrey+($halfnpntsy*$2))) $(($centrez+($halfnpntsz*$2))) )
corners+=( $(($centrex-($halfnpntsx*$2))) $(($centrey+($halfnpntsy*$2))) $(($centrez+($halfnpntsz*$2))) )

p=0
for ((n=0; n < $((${#corners[*]}/3)); n++)); do
  printf "ATOM      $((n+1))  CO$n BOX   1      " >> ${1%.*}_box.pdb
  for ((c=0; c < 3; c++)); do
    print -f  "%8.3f" -- "${corners[$p]}" >> ${1%.*}_box.pdb
    ((p++))
  done
  print >> ${1%.*}_box.pdb
done
print "CONECT    1    2    4    5
CONECT    2    1    3    6
CONECT    3    2    4    7
CONECT    4    1    3    8
CONECT    5    1    6    8
CONECT    6    2    5    7
CONECT    7    3    6    8
CONECT    8    4    5    7" >> ${1%.*}_box.pdb
}



function makegpf {
if [[ ${#arglist[*]} != 3 ]]; then
  error="USAGE: autoAD.ksh -g ligand.pdbqt protein_H.pdbqt pocket.[mol2/pdb]"
  errhan
fi
if [[ ${#arglist[*]} > 2 ]] && [[ ${arglist[1]##*.} == pdb ]]; then
  argno=1
  if [[ -e ${arglist[0]%.*}1.pdbqt ]]; then
    ligand=${arglist[0]%.*}1.pdbqt
  else
    ligand=${arglist[0]%.*}.pdbqt
  fi
  input_receptor=${arglist[${argno}]%.*}_H.pdbqt
elif [[ ${#arglist[*]} > 2 ]] && [[ ${arglist[1]##*.} == pdbqt ]] &&  [[ ${arglist[0]##*.} == pdbqt ]]; then
  ligand=${arglist[0]%.*}.pdbqt
  input_receptor=${arglist[1]%.*}.pdbqt
  argno=1
elif [[ ${arglist[0]##*.} == pdb?(qt) ]]; then
  if [[ ${arglist[1]##*.} != pdb ]]; then
      print "Converting ${arglist[1]} to ${arglist[1]%.*}.pdb"
      /usr/local/bin/babel -i${arglist[1]##*.} ${arglist[1]} -opdb ${arglist[1]%.*}.pdb
  fi
  argno=0
  ligand=${arglist[1]%.*}.pdb
  input_receptor=${arglist[${argno}]%%?(_H).*}_H.pdbqt
elif [[ ${arglist[0]##*.} == pdbqt ]]; then
  argno=0
  input_receptor=${arglist[${argno}]}
elif [[ ${arglist[1]##*.} == pdbqt ]]; then
  argno=1
  input_receptor=${arglist[1]}
  if [[ -e ${arglist[0]%.*}1.pdbqt ]]; then
    ligand=${arglist[0]%.*}1.pdbqt
  else
    ligand=${arglist[0]%.*}.pdbqt
  fi
else
  error="Are you sure your arguments are in the right order? autoAD.ksh -g ligand.pdbqt protein_H.pdbqt pocket.[mol2/pdb]"
  errhan
fi
rm -f *~
print "Making grid parameter file for ${arglist[${argno}]%.*}"
print "Ligand is ${arglist[$((argno+1))]}"
print "parameter_file /usr/people/douglas/programs/autodock/AD4.1_bound.dat
$(boxsize ${arglist[$((argno+1))]} $spacing)
spacing $spacing" > gpftemplate_tmp 
print "Box written to ${arglist[$((argno+1))]%.*}_box.pdb"
print "~douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l $ligand -r $input_receptor -i gpftemplate_tmp -d ./ "
/usr/people/douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -l $ligand -r $input_receptor -i gpftemplate_tmp -d ./ 
mv -f ${arglist[${argno}]%%?(_H).*}*.gpf ${arglist[${argno}]%%?(_H).*}_H.gpf
sed -i "s/npts.*/$(grep npts gpftemplate_tmp)                        # num.grid points in xyz/" ${arglist[${argno}]%%?(_H).*}_H.gpf
atomtypes=$(grep "^map .*\.map" ${arglist[${argno}]%%?(_H).*}_H.gpf | wc | awk '{print $1}')
if (( atomtypes > 14 )); then
  print "Too many atom types, splitting gpf files ..."
  cp ${arglist[${argno}]%%?(_H).*}_H.gpf ${arglist[${argno}]%%?(_H).*}_H2.gpf
  maps=($(grep "^map .*\.map" ${arglist[${argno}]%%?(_H).*}_H.gpf | awk '{print $2}'))
  for ((n=0; n < $((atomtypes-14)); n++)); do
    sed -i "/.*${maps[$n]}.*/d" ${arglist[${argno}]%%?(_H).*}_H.gpf
    deltype=${maps[$n]%.map}; deltype=${deltype#*.}
    deltypes+="$deltype "
  done
  for deltype in $deltypes; do
    sed -i "s/\(^ligand_types \)$deltype \(.*\)/\1\2/" ${arglist[${argno}]%%?(_H).*}_H.gpf 
  done
  for ((n=$n; n < $atomtypes; n++)); do
    sed -i "/${maps[$n]}/d" ${arglist[${argno}]%%?(_H).*}_H2.gpf
  done
  sed -i  "s/^ligand_types.*/ligand_types $deltypes # ligand atom types/" ${arglist[${argno}]%%?(_H).*}_H2.gpf
fi
#rm -f gpftemplate_tmp
}


function makemaps {
if [[ ${#arglist[*]} > 1 ]] && [[ ${arglist[1]##*.} == pdb?(qt) ]]; then
  argno=1
elif [[ ${arglist[0]##*.} == pdb ]]; then
  argno=0
else
  error="Are you sure your arguments are in the right order?"
  errhan
fi
print "Calculating grids, logfile: ${arglist[${argno}]%%?(_H).*}_H.glg"
if [[ -e ${arglist[${argno}]%%?(_H).*}_H2.gpf ]]; then
  print "Split grid parameter files detected, running Autogrid twice ..."
  /usr/people/douglas/programs/autodock/autogrid4 -p ${arglist[${argno}]%%?(_H).*}_H2.gpf | tee ${arglist[${argno}]%%?(_H).*}_H2.glg
fi
sleep 3
print "/usr/people/douglas/programs/autodock/autogrid4 -p ${arglist[${argno}]%%?(_H).*}_H.gpf | tee ${arglist[${argno}]%%?(_H).*}_H.glg"
sleep 3
/usr/people/douglas/programs/autodock/autogrid4 -p ${arglist[${argno}]%%?(_H).*}_H.gpf  2>&1| tee ${arglist[${argno}]%%?(_H).*}_H.glg
}


function makedpfs {
if [[ -n $max_params ]]; then
  ga_runs=20
  function calculate_minimizations { 
  print "Reading number of torsions from $file"; sleep 3
  torsions=$(tail -1 $file | cut -c8-11 | tr -d "[:space:]")
  if (( torsions < 11 )); then
    ((ga_num_evals=(((${torsions}**2)*0.2627)+(0.1551*${torsions})+0.2827)*1000000 ))
  else
    ((ga_num_evals=((${torsions}*1.1364)+13.636)*1000000 ))
  fi
  }
else
  ga_runs=10
  minimizations=(125000 250000 500000 1000000 2000000 4000000 5500000 7000000 8000000 9000000 10000000)
  function calculate_minimizations {
  ga_num_evals=${minimizations[$(head -1 $file | cut -c7-10 | tr -d "[:space:]")]-10000000}
  }
fi


if [[ ${#arglist[*]} < 2 ]]; then
  error="Specify compoundlist.sdf AND receptor.pdb"
  errhan
fi
print "parameter_file /usr/people/douglas/programs/autodock/AD4.1_bound.dat 
ga_pop_size 300
ga_num_generations 30000
ga_run $ga_runs" > dpftemplate_tmp


for file in $(ls -1 ${ligid}*.pdbqt | sort -n --key=1.$sortkey); do
  calculate_minimizations
  print "Making docking parameter file for ${file}"
  if (( pose == 0 )); then
    print "/usr/people/douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -i dpftemplate_tmp -l $file -r ${arglist[1]%%?(_H).*}_H.pdbqt -p ga_num_evals=$ga_num_evals"
    /usr/people/douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -i dpftemplate_tmp -l $file -r ${arglist[1]%%?(_H).*}_H.pdbqt -p ga_num_evals=$ga_num_evals > /dev/null
  else
    print "/usr/people/douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -i dpftemplate_tmp -l $file -r ${arglist[1]%%?(_H).*}_H.pdbqt -p ga_num_evals=$ga_num_evals -s"
    /usr/people/douglas/programs/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf4.py -i dpftemplate_tmp -l $file -r ${arglist[1]%%?(_H).*}_H.pdbqt -p ga_num_evals=$ga_num_evals -s > /dev/null
    sed -i "s/^ga_run 10/ga_run 2/" ${file%.pdb?(qt)}_${arglist[1]%%?(_H).*}_H.dpf
  fi
done

sleep 5
sed -i "s/^unbound_model extended/unbound_model bound/" *.dpf
sed -i "s/# / # /" *.dpf
sed -i "s/outlev 1/outlev 2/" *.dpf

#rm -f dpftemplate_tmp
}



function prepDOCK6 {
print "DOCK mode activated"
jobaddon=_DOCKed
if [[ ! -e hostnames.txt ]]; then
  print "$hostnames" > hostnames.txt
fi

rm -f OUTSPH temp* protein.sph

#  argno=1
#function preprep {
if [[ ${#arglist[*]} > 1 ]] && [[ ${arglist[1]##*.} == pdb ]]; then
  print argno is 1
  argno=1
elif [[ ${arglist[0]##*.} == pdb ]]; then
  print argno is 0
  argno=0
else
  error="Are you sure your arguments are in the right order? Should be: library.sdf protein.pdb bindsite.mol2"
  errhan
fi


print "Preparing receptor ... "
grep -v "^TER" ${arglist[${argno}]} | grep -v "^REMARK" > ${arglist[${argno}]%.*}_tmp.$$
OLDPYTHONHOME=$PYTHONHOME
export PYTHONHOME=/usr/people/douglas/programs/Python-2.7.11/build/

print "Hydrogenating and ionizing ..."
print "/usr/people/douglas/programs/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py --ff=PARSE --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr"
/usr/people/douglas/programs/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py --ff=PARSE --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr

#unset PYTHONHOME
#print "usr/people/douglas/programs/pdb2pqr-1.8/pdb2pqr.py --ff=PARSE --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr"
#/usr/people/douglas/programs/pdb2pqr-1.8/pdb2pqr.py --ff=PARSE --with-ph=7.4  ${arglist[${argno}]%.*}_tmp.$$ ${arglist[${argno}]%.*}_H.pqr
export PYTHONHOME=$OLDPYTHONHOME

if [[ ! -s ${arglist[${argno}]%.*}_H.pqr ]]; then
  print "WARNING: PDB2PQR has failed, switching to OpenBabel..."
  print /usr/local/bin/babel -ipdb ${arglist[${argno}]%.*}_tmp.$$ -opdb ${arglist[${argno}]%.*}_H.pqr -p 7.4
  /usr/local/bin/babel -ipdb ${arglist[${argno}]%.*}_tmp.$$ -opdb ${arglist[${argno}]%.*}_H.pqr -p 7.4
fi
rm -f ${arglist[${argno}]%.*}_tmp.$$

/usr/local/bin/babel -ipqr ${arglist[${argno}]%.*}_H.pqr -omol2 ${arglist[${argno}]%.*}_H.mol2

print "@<TRIPOS>SUBSTRUCTURE
  1      GARBAGE 1       PERM    0       ****    ****    0       ROOT" >> ${arglist[${argno}]%.*}_H.mol2
edits=($(head -3 ${arglist[${argno}]%.*}_H.mol2 | tail -1))
sed -i "s/${edits[*]}/ ${edits[0]} ${edits[1]} 1 ${edits[3]}/" ${arglist[${argno}]%.*}_H.mol2

/usr/local/bin/babel -ipqr ${arglist[${argno}]%.*}_H.pqr -opdb ${arglist[${argno}]%.*}_noH.pdb -d
print "Checking library SD file for errors ..."
dot=$(head -1 ${arglist[0]} | cut -f1 -d" ")
sed -i "s/^$dot/${dot#.}/" ${arglist[0]} 
sed -i "s/0999 V2000/0   1 V2000/" ${arglist[0]} 
print Removing any hydrogen atoms from SD file ...
/usr/local/bin/babel -isdf ${arglist[0]} -osdf ${ligid}_noH_tmp.sdf.$$ -d
print Adding hydrogen atoms to compounds and converting to mol2 ...
print "/usr/local/bin/babel -m -isdf ${ligid}_noH_tmp.sdf.$$ -omol2 ${ligid}.mol2 -p 7.4"
/usr/local/bin/babel -m -isdf ${ligid}_noH_tmp.sdf.$$ -omol2 ${ligid}.mol2 -p 7.4
rm -f ${ligid}_noH_tmp.sdf.*


#Get only the protein residues that are close to the ligand to make dms and sphgen job smaller
/usr/local/bin/babel -imol2 ${arglist[2]} -opdb ${arglist[2]%.*}.pdb

print "~douglas/scripts/bindsite.ksh ${arglist[${argno}]%.*}.pdb ${arglist[2]%.*}.pdb 10"
~douglas/scripts/bindsite.ksh ${arglist[${argno}]%.*}.pdb ${arglist[2]%.*}.pdb 10

print "~douglas/programs/DMS/bin/dms ${arglist[${argno}]%.*}_site.pdb -n -w1.4 -v -o ${arglist[${argno}]%.*}_noH.ms"
~douglas/programs/DMS/bin/dms ${arglist[${argno}]%.*}_site.pdb -n -w1.4 -v -o ${arglist[${argno}]%.*}_noH.ms

rm -f INSPH OUTSPH temp* ${arglist[${argno}]%.*}_noH_site.sph
print "${arglist[${argno}]%.*}_noH.ms
R
X
0.0
4.0
1.4
${arglist[${argno}]%.*}_noH_site.sph" > INSPH; wait


~douglas/programs/dock6/bin/sphgen 

~douglas/programs/dock6/bin/sphere_selector ${arglist[${argno}]%.*}_noH_site.sph ${arglist[2]} 5
mv selected_spheres.sph ${arglist[${argno}]%.*}_selected_spheres.sph

print "Y
10
${arglist[${argno}]%.*}_selected_spheres.sph
1
${arglist[${argno}]%.*}_noH_box.pdb" > box.in; wait

~douglas/programs/dock6/bin/showbox < box.in

print "${arglist[${argno}]%.*}_selected_spheres.sph
1
N
sphgen_cluster.pdb" > showsphere.in

~douglas/programs/dock6/bin/showsphere < showsphere.in

print "compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  yes
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  ${arglist[${argno}]%.*}_H.mol2
box_file                       ${arglist[${argno}]%.*}_noH_box.pdb
vdw_definition_file            /usr/people/douglas/programs/dock6/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              grid
contact_cutoff_distance        4.5" > grid.in; wait

if [[ $(head -5 ${arglist[${argno}]%.*}_H.mol2 | tail -1) == GASTEIGER ]]; then
  rm -f tmp
  print "$(head -4 ${arglist[${argno}]%.*}_H.mol2)\n$(tail --lines=+6 ${arglist[${argno}]%.*}_H.mol2)" >> tmp
  mv tmp ${arglist[${argno}]%.*}_H.mol2
fi

~douglas/programs/dock6/bin/grid -i grid.in
}



function makeDOCK6dpf {
if [[ ${#arglist[*]} < 2 ]]; then
  error="Specify compoundlist.sdf AND receptor.pdb"
  errhan
fi

for file in $(ls -1 ${ligid}+([[:digit:]]).mol2 | sort -n --key=1.$sortkey); do
  print "Making docking parameter file for ${file}"
  print "conformer_search_type                                        flex
ligand_atom_file                                             $file
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ${arglist[1]%%?(_H).*}_selected_spheres.sph
max_orientations                                             20000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
flexible_ligand                                              yes
user_specified_anchor                                        no 
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            yes
clash_overlap                                                0.5
write_growth_tree                                            no
bump_filter                                                  yes
bump_grid_prefix                                             grid
max_bumps_anchor                                             12
max_bumps_growth                                             12
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       grid
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_tors_premin_iterations                               10
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                1000
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                /usr/people/douglas/programs/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /usr/people/douglas/programs/dock6/parameters/flex.defn
flex_drive_file                                              /usr/people/douglas/programs/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${file%.mol2}_DOCKed
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
multigrid_score_secondary                                    no
SASA_score_secondary                                         no" > ${file%.mol2}.dockin

done
}


function autodock {
if [[ ${#arglist[*]} > 1 ]]; then
  arg=${arglist[0]%.*}1_${arglist[1]%.pdb}_H.dpf
else
  arg=${arglist[0]}
fi
#print "Use tail -f ${arg%.dpf}.dlg to monitor"
print "Running \"/usr/people/douglas/programs/autodock/autodock4 -p $arg -l ${arg/%.dpf/.dlg} > /dev/null 2>&1 &;\""
trap 'print "Killed autodock"; kill $pid' INT 
/usr/people/douglas/programs/autodock/autodock4 -p $arg -l ${arg/%.dpf/.dlg} > /dev/null 2>&1 &;
pid=$(jobs -p %+)
sleep 3
tail -f --pid=$pid ${arg%.dpf}.dlg

}


function listvina {
if [[ ! -e vina.dpf ]]; then
  error="Vina.dpf not found; please run autoAD.ksh -j to generate it."
  errhan $0
fi
centers=($(grep gridcenter vina.dpf))
npts=($(grep npts vina.dpf))
receptor=$(grep .pdbqt vina.dpf)
print $receptor
for file in $(< ${arglist[0]}); do
  print $file
  if ! grep -q -s "Writing output ... done." ${file%.dpf}_out.dlg; then
    print Running Vina ...
    /localdisk/software/gwovina-1.0/bin/gwovina --receptor $receptor --ligand $file --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) $max_params | tee -a ${arglist[0]%.pdbqt}_out.dlg
    print run_xscore ${receptor} ${arglist[0]%.pdbqt}_out.pdbqt
    run_xscore ${receptor} ${arglist[0]%.pdbqt}_out.pdbqt
  else
    print "${file%.dpf}.dlg already complete!"
  fi
done
}


function listdock {
for file in $(< ${arglist[0]}); do
  if ! grep -q -s "Successful Completion" ${file%.dpf}.dlg; then
    /usr/people/douglas/programs/autodock/autodock4 -p ${file%.dpf}.dpf > ${file%.dpf}.dlg 2>&1 
    print "Finished $file"
  else
    print "${file%.dpf}.dlg already complete!"
  fi
done
}


function dyndock {
dynend=0
while (( dynend < 1 )); do
  if [[ $(grep "^Finished" ${arglist[0]} 2> /dev/null) != Finished ]] && [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) != Die ]]; then
    for file in $(< ${arglist[0]}); do
      if [[ $(grep "Successful Completion" ${file%.dpf}.dlg 2> /dev/null) != "Successful Completion" ]] && [[ {$file%.dpf} != Finished ]]; then
        print "Dyndock starting ${file%.dpf}.dpf on $(hostname -s); dynend is $dynend; date is $(date); jobs read from ${arglist[0]}" >> ${arglist[0]%.txt}.log
        /usr/people/douglas/programs/autodock/autodock4 -p ${file%.dpf}.dpf > ${file%.dpf}.dlg; wait
      fi
    done
    until [[ $(cat ${arglist[0]} 2> /dev/null) == Finished ]]; do
      print "Finished" > ${arglist[0]}
    done
  fi
  sleep 1
  if [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) == Die ]]; then
    print "Die command word detected; contents of ${arglist[0]} is $(cat ${arglist[0]}), dynend is $dynend" >> ${arglist[0]%.txt}.log
    ((dynend++))
  fi
done
print "Dyndock finished ${arglist[0]} on $(hostname -s); dynend is $dynend; date is $(date); last ligand was ${file%.dpf}" >> ${arglist[0]%.txt}.log
}

function vinadpf { 
print "Making Vina parameter file vina.dpf"
if (( ${#arglist[*]} < 2 )); then
  error="Requires 2 arguments: receptor_withHs.pdbqt, binding_site_definition.[mol2|pdb]."
  errhan $0
else
  print "boxsize ${arglist[${#arglist[*]}-1]} $spacing > vina.dpf"
  boxsize ${arglist[${#arglist[*]}-1]} $spacing > vina.dpf
  print "Box written to ${arglist[${#arglist[*]}-1]%.*}_box.pdb"
  print ${arglist[${#arglist[*]}-2]/%.pdb/_H.pdbqt} >> vina.dpf
fi
print $ligid >> vina.dpf
print ${arglist[${#arglist[*]}-1]} >> vina.dpf
print -- "$max_params" >> vina.dpf
print -- "$num_modes" >> vina.dpf
}

function vina_oneoff {
if [[ ! -e vina.dpf ]]; then
  error="Vina.dpf not found; please run autoAD.ksh -j to generate it."
  errhan $0
fi
if (( ${#arglist[*]} != 1 )); then
  error="-n takes 1 argument - the ligand.pdbqt - the rest it gets from vina.dpf"
  errhan
fi
centers=($(grep gridcenter vina.dpf))
npts=($(grep npts vina.dpf))
receptor=$(grep .pdbqt vina.dpf)
max_params=$(grep exhaustiveness vina.dpf)
num_modes=$(grep num_modes vina.dpf)
if (( pose == 1 )); then
  prepose="--local_only"
fi 
rm -f ${arglist[0]%.pdbqt}_out.dlg
print Running Vina ...
print "/localdisk/software/gwovina-1.0/bin/gwovina --receptor $receptor --ligand ${arglist[0]} --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) $max_params $num_modes $prepose | tee -a ${arglist[0]%.pdbqt}_out.dlg"
/localdisk/software/gwovina-1.0/bin/gwovina --receptor $receptor --ligand ${arglist[0]} --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) $max_params $num_modes $prepose | tee -a ${arglist[0]%.pdbqt}_out.dlg

#print run_xscore ${receptor} ${arglist[0]%.pdbqt}_out.pdbqt
#run_xscore ${receptor} ${arglist[0]%.pdbqt}_out.pdbqt
}


function dynvina {
dynend=0
while (( dynend < 1 )); do
  if [[ $(grep "^Finished" ${arglist[0]} 2> /dev/null) != Finished ]] && [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) != Die ]]; then
    for ligand in $(< ${arglist[0]}); do
      if [[ $(grep "Writing output ... done." ${ligand%.pdbqt}.dlg 2> /dev/null) != "Writing output ... done." ]]; then
        print "Dyndock starting $ligand on $(hostname -s); dynend is $dynend; date is $(date); jobs read from ${arglist[0]}" >> ${arglist[0]%.txt}.log
        centers=($(grep gridcenter vina.dpf)) 
        npts=($(grep npts vina.dpf))
        receptor=$(grep _H.pdbqt vina.dpf)
        max_params=$(grep exhaustiveness vina.dpf)
        print Running Vina ...
        print "Vina docking $ligand using ${arglist[0]} on $(hostname -s) max_params is $max_params" > ${ligand%.pdbqt}_out.dlg
        print "/localdisk/software/gwovina-1.0/bin/gwovina --receptor $receptor --ligand $ligand --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) $max_params" >> ${ligand%.pdbqt}_out.dlg
        /localdisk/software/gwovina-1.0/bin/gwovina --receptor $receptor --ligand $ligand --center_x ${centers[1]} --center_y ${centers[2]} --center_z ${centers[3]} --size_x $(((${npts[1]}+1)*$spacing)) --size_y $(((${npts[2]}+1)*$spacing)) --size_z $(((${npts[3]}+1)*$spacing)) $max_params >> ${ligand%.pdbqt}_out.dlg 2>&1 ; wait
        checks=0
        until [[ $(grep "Writing output ... done." ${ligand%.pdbqt}_out.dlg) == "Writing output ... done." || $(grep "^Parse error" ${ligand%.pdbqt}_out.dlg) != "" ]]; do
          sleep 1
          ((checks++))
          if (( checks > 240 )); then
            failed=1
            break
          fi
        done 
        if (( failed == 1 )); then
          rm -rf ${ligand%.pdbqt}_out.dlg
        else
          print "Vina done docking $ligand using ${arglist[0]} on $(hostname -s)" >> ${ligand%.pdbqt}_out.dlg
        fi
        failed=0
#Run X-score on output pose
        run_xscore ${receptor} ${ligand%.pdbqt}_out.pdbqt
#        run_drugscore ${receptor} ${ligand%.pdbqt}_out.pdbqt
      fi
    done
    until [[ $(cat ${arglist[0]} 2> /dev/null) == Finished ]]; do
      print "Finished" > ${arglist[0]}
    done
  fi
  sleep 1
  if [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) == Die ]]; then
    print "Die command word detected; contents of ${arglist[0]} is $(cat ${arglist[0]})" >> ${arglist[0]%.txt}.log
    ((dynend++))
  fi
done
print "Dynvina finished ${arglist[0]} on $(hostname -s); dynend is $dynend; date is $(date); last ligand was $ligand" >> ${arglist[0]%.txt}.log
}


function dynDOCK6 {

dynend=0
while (( dynend < 1 )); do
  if [[ $(grep "^Finished" ${arglist[0]} 2> /dev/null) != Finished ]] && [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) != Die ]]; then
    for file in $(< ${arglist[0]}); do
      if [[ $(grep "Total elapsed time" ${file%.dpf}_DOCKed.dlg 2> /dev/null) != "Total elapsed time" ]] && [[ {$file%.dpf} != Finished ]]; then
        print "dynDOCK6 starting ${file%.dockin}.dockin on $(hostname -s); dynend is $dynend; time is $(date); jobs read from ${arglist[0]}" >> ${arglist[0]%.txt}.log
        ~douglas/programs/dock6/bin/dock6 -i ${file%.dockin}.dockin > ${file%.dockin}_DOCKed.dlg; wait
      fi
    done
    until [[ $(cat ${arglist[0]} 2> /dev/null) == Finished ]]; do
      print "Finished" > ${arglist[0]}
    done
  fi
  sleep 1
  if [[ $(grep "^Die" ${arglist[0]} 2> /dev/null) == Die ]]; then
    print "Die command word detected; contents of ${arglist[0]} is $(cat ${arglist[0]}), dynend is $dynend" >> ${arglist[0]%.txt}.log
    ((dynend++))
  fi
done
print "dynDOCK6 finished ${arglist[0]} on $(hostname -s); dynend is $dynend; date is $(date); last ligand was ${file%.dpf}" >> ${arglist[0]%.txt}.log
}



function autodist {
#.Script to distribute Autodock jobs to more than one cluster PC 
batch=1

case $mode_selection in 
  y) ligid=$(cat ligid.txt)
     print Autodist running Autodock - ligID: $ligid
     jobs=($(ls *.dpf))
     if [[ -n $(ls *.dlg 2> /dev/null) ]]; then
       print "Previous docking logfiles detected, please wait ..."
       for done in $(tail *.dlg | grep -B6 Success | grep .dlg | cut -f2 -d" "); do
#        jobs=(${jobs[*]/${done/%.dlg/.dpf}}); unset jobs[${#jobs[*]}-1]
         oldjobs=(${jobs[*]})
         jobs=(${oldjobs[*]//${done/%.dlg/.dpf}})
       done       
     fi
     jobsno=${#jobs[*]}
     extension=txt
     if [[ ! -e hostnames.txt ]] || [[ -z $jobs ]]; then
       error="Cannot find input files"
       errhan
     fi ;;
  z) ligid=$(head -4 vina.dpf | tail -1); ligid=${ligid%.sdf}
     print Autodist running Vina - ligID: $ligid
     jobaddon=_out
     jobs=($(ls -1 ${ligid}!(*_out).pdbqt | sort -n --key=1.$sortkey))
     if [[ -n $(ls *.dlg 2> /dev/null) ]]; then
       print "Previous docking logfiles detected, please wait ..."
       for done in $(grep "Writing output ... done." *.dlg | cut -f1 -d":"); do
#        jobs=(${jobs[*]//${done/%_out.dlg/.pdbqt}}); unset jobs[${#jobs[*]}-1]
         oldjobs=(${jobs[*]})
         jobs=(${oldjobs[*]//${done/%_out.dlg/.pdbqt}})
       done       
     fi
     jobsno=${#jobs[*]}
     extension=txt
     if [[ ! -e hostnames.txt ]] || [[ -z $jobs ]]; then
       error="Cannot find input files"
       errhan
     fi ;;
  f) print Autodist running DOCK 6.8 - ligID: $ligid
     jobaddon=_DOCKed
     jobs=($(ls -1 ${ligid}[0-9]*.dockin | sort -n --key=1.$sortkey))
     if [[ -n $(ls *.dlg 2> /dev/null) ]]; then
       print "Previous docking logfiles detected, please wait ..."
       for done in $(grep "Total elapsed time" *.dlg | cut -f1 -d":"); do
         oldjobs=(${jobs[*]})
         jobs=(${oldjobs[*]//${done/%${jobaddon}.dlg/.dockin}})
       done       
     fi
     jobsno=${#jobs[*]}
     print $jobsno jobs left to do
     extension=txt
     if [[ ! -e hostnames.txt ]] || [[ -z $jobsno ]]; then
       error="Cannot find input files"
       errhan
     fi  ;;
esac


#Get all docking jobs
nodeinfo=(nodes=($(cut -f1 hostnames.txt)) cpuscores=($(cut -f2 hostnames.txt)))
print "Nodes: ${nodeinfo.nodes[*]}"


function dynamic {


totaljobs=${#jobs[*]}
print "$totaljobs jobs will be farmed out."
alljobs=(${jobs[*]})
rm -f docklist_*.log

while [[ -n $jobs ]]; do
  unset end
  unset lastnode
  nodesdone=0
  rm -f docklist_*_*.txt docklist_*_*.mol2 autodist_*_tmp; wait
  while [[ -z $end ]]; do
#Assign jobs to each node according to batch size
    jobs_assigned=0
    c=0
    for node in ${nodeinfo.nodes[*]}; do
      ((c++))
      for ((n=0; n < batch; n++)); do
        if ((jobs_assigned < totaljobs)); then
          print ${jobs[$jobs_assigned]/%.dpf} >> docklist_${node}_${c}.$extension
          ((jobs_assigned++))
        fi
      done
    done

#Clean up node scripts
    for node in ${nodeinfo.nodes[*]}; do
      rm -f autodist_${node}_tmp
    done
    sleep 3
#Login to each node and execute autoAD.ksh using the docking list as input
    c=0
    for node in ${nodeinfo.nodes[*]}; do
      ((c++))
      print "/localdisk/home/s2173344/Research_Project/autoAD.ksh -${mode_selection} docklist_${node}_${c}.$extension&" >> autodist_${node}_tmp
      if [[ ! -e autodist_${node}_tmp ]]; then
        sleep 5
      fi
      chmod +x autodist_${node}_tmp
    done
    
    startedjobs=0
    for node in ${nodeinfo.nodes[*]}; do
      if [[ $node != $lastnode ]] && (( startedjobs < totaljobs )); then
        print Logging in to $node to run: "cd $(pwd); nice +9 nohup ./autodist_${node}_tmp  </dev/null >&/dev/null &; exit"
        ((startedjobs++))
        priority=9
        if [[ $(whoami) != douglas ]] || [[ $node == itiwks && ${mode_selection} == z ]]; then
          priority=18
        fi
        ./autodist_${node}_tmp >autodist.log 2>&1 </dev/null &
      fi
      lastnode=$node
    done

#Now keep checking the docklist files to monitor each node's progress
    n=0
    sleep 5
    print "Dir: $(pwd)"
    if [[ -z $verbose ]]; then
      printf "\nAutodist running [ "
    else
      print "\nAutodist running in verbose mode ..." 
    fi
    sleep 5
    while (( end < 1 )); do
      for ((l=0; l < 4; l++))
      do
        sleep 1
        if [[ -z $verbose ]]; then
          printf "."
        fi
      done
      if [[ -z $verbose ]]; then
        printf '\010\010\010\010    ]\010\010\010\010\010'
      fi
      c=0
      for node in ${nodeinfo.nodes[*]}; do
        ((c++))
        if [[ $(grep "^Finished" docklist_${node}_${c}.$extension 2> /dev/null) == Finished ]]; then
          if ((jobs_assigned < totaljobs)) || (( $(wc -l  docklist_${node}_${c}.$extension | awk '{print $1}') != 1 )); then
            until [[ ! -e  docklist_${node}_${c}.$extension ]]; do
              rm -fr docklist_${node}_${c}.$extension
            done
          fi
#What can happen now is the node script prints Finished to the file that has just been deleted
          for ((n=0; n < batch; n++)); do
            if ((jobs_assigned < totaljobs)); then
#              if [[ ! -e docklist_${node}_${c}.$extension ]]; then
                print ${jobs[$jobs_assigned]/%.dpf} >> docklist_${node}_${c}.$extension
                ((jobs_assigned++))
#              fi
            fi
          done

        fi
      done

#Check to see there's no jobs left and that each node has finished, if so stop autoAD on that node
      c=0
      if ((jobs_assigned >= totaljobs)); then
        for node in ${nodeinfo.nodes[*]}; do
          ((c++))
          if [[ $(cat docklist_${node}_${c}.$extension  2> /dev/null) == Finished ]]; then
            if [[ -n $verbose ]]; then
              print "Printing Die to docklist_${node}_${c}.$extension since jobs_assigned ($jobs_assigned) is >= totaljobs ($totaljobs) and keyword $(cat docklist_${node}_${c}.$extension) was found."
            fi
            write_tries=1 
            until [[ $(cat docklist_${node}_${c}.$extension) == Die ]] || (( write_tries > 100 )); do
              print "Die" > docklist_${node}_${c}.$extension
              ((write_tries++))
            done
            ((nodesdone++))
          fi
        done
        c=0
        if [[ -n $verbose ]]; then
          print "Number of docklists is $(ls docklist_*_*.txt | wc -w), nodesdone is $nodesdone."
        fi
        if ((nodesdone >= $(ls docklist_*_*.txt | wc -w))); then
          ((end++))
        fi
      fi
    done
  done
  rm -f Finished${jobaddon}.dlg Die${jobaddon}.dlg
  sleep 60
  unset jobs
  for job in ${alljobs[*]}; do
    if [[ ! -e ${job%.*}${jobaddon}.dlg ]]; then
      jobs+=(${job})
    fi
  done
  printf "Jobs found: $totaljobs  Jobs attempted: $jobs_assigned  Jobs failed: ${#jobs[*]} ]\n" 
  if (( ${#jobs[*]} != 0 )); then
    print "Will attempt the following jobs again: ${jobs[*]}"
  fi
  totaljobs=${#jobs[*]}
done
rm -f Finished${jobaddon}.dlg Die${jobaddon}.dlg
}


function oneoff {
print ${nodeinfo.cpuscores[*]}

if [[ $(whoami) == douglas ]]; then
  priority=9
else
  priority=18
fi

for score in ${nodeinfo.cpuscores[*]}; do
  ((total_cpu_score+=score))
done
print Total cpu score is $total_cpu_score

((score_multiplier=1.0/${total_cpu_score}))
print Score multiplier is $score_multiplier

typeset -A node_multipliers
for node in ${nodeinfo.nodes[*]}; do
  ((portion=${nodeinfo.cpuscores[$n]}*score_multiplier))
  node_multipliers+=( [$node]=$portion )
  ((n++))
done
n=

typeset -A adjusted_batchsizes
for node in ${nodeinfo.nodes[*]}; do
  print Adjusted batchsize for $node is $((${jobsno}*${node_multipliers[$node]}))
  adjusted_batchsizes+=( [$node]=$((${jobsno}*${node_multipliers[$node]})) )
  ((total+=${adjusted_batchsizes[$node]%.*}))
  ((sum+=${node_multipliers[$node]}))
done
print $sum should be close to 1
((remainder=${jobsno}-total))

#Divvy up the jobs to each node, then create a list of dockings to perform for each node
print Number of jobs is ${jobsno}
print Number of nodes is ${#nodeinfo.nodes[*]}
print Remainder is $remainder

c=
n=
rm -f docklist_*_*.*
  start=1
  for node in ${nodeinfo.nodes[*]}; do
    ((c++))
    if [[ $jobaddon! = _DOCKed ]]; then
      while (( n < ${adjusted_batchsizes[$node]%.*} )); do
        print ${jobs[$o]/%.dpf} >> docklist_${node}_${c}.$extension
        ((n++))
        ((o++))
      done
      if (( remainder > 0 )); then
        print ${jobs[$o]/%.dpf} >> docklist_${node}_${c}.$extension
        ((remainder--))
        ((o++))
        p=1
      else
        p=0
      fi
    else
      max=$(grep -n '\@<TRIPOS>MOLECULE' ${ligid}_OBH.mol2 | wc -l)
      ((end+=${adjusted_batchsizes[$node]%.*}))
      startline_no=$(grep -n '\@<TRIPOS>MOLECULE' ${ligid}_OBH.mol2 | head -$start | tail -1 | cut -f1 -d":")
      endline_no=$(grep -n '\@<TRIPOS>MOLECULE' ${ligid}_OBH.mol2 | head -$((end+1)) | tail -1 | cut -f1 -d":"); ((endline_no--))
      if (( end != max )); then
        ((tailinput=endline_no-startline_no+1))
      else
        endline_no=$(cat ${ligid}_OBH.mol2 | wc -l)
        ((tailinput=endline_no-startline_no+1))
      fi
#     print "head -$endline_no ${ligid}_OBH.mol2 | tail -$tailinput $endline_no $startline_no"
      head -$endline_no ${ligid}_OBH.mol2 | tail -$tailinput > docklist_${node}_${c}.mol2
      ((n=end-start+1))
      ((start+=${adjusted_batchsizes[$node]%.*}))
      if (( remainder > 0 )); then
        startline_no=$(grep -n '\@<TRIPOS>MOLECULE' ${ligid}_OBH.mol2 | head -$start | tail -1 | cut -f1 -d":")
        endline_no=$(grep -n '\@<TRIPOS>MOLECULE' ${ligid}_OBH.mol2 | head -$((start+1)) | tail -1 | cut -f1 -d":")
        head -$((endline_no-1)) ${ligid}_OBH.mol2 | tail -$((endline_no-startline_no+2)) >> docklist_${node}_${c}.mol2
        ((start++))
        ((end++))
        ((remainder--))
        p=1
      else
        p=0
      fi
      ((written+=n+p))
      print Written $((n+p)) jobs to docklist_${node}_${c}; 
      n=0
    fi
  done
  print ${jobsno} jobs were found, $written have been assigned
  c=0
  

#Clean up node scripts
for node in ${nodeinfo.nodes[*]}; do
    rm -f autodist_${node}_tmp
done

#Login to each node and execute autoAD.ksh using the docking list as input
for node in ${nodeinfo.nodes[*]}; do
  ((c++))
  print "/usr/people/douglas/scripts/docking/autoAD.ksh -${mode_selection} ./docklist_${node}_${c}.$extension&" >> autodist_${node}_tmp
  until [[ -e autodist_${node}_tmp ]]; do
    sleep 1
  done
  chmod +x autodist_${node}_tmp
done

for node in ${nodeinfo.nodes[*]}; do
  if [[ $node != $lastnode ]]; then
    print Logging in to $node to run: 'ssh $node "cd $(pwd); nice +${priority} nohup ./autodist_${node}_tmp  </dev/null >&/dev/null &; exit"'
    ssh $node "cd $(pwd); nice +${priority} nohup ./autodist_${node}_tmp  </dev/null >&/dev/null &; exit"
  fi
  lastnode=$node
done
}


if [[ $1 == "-d" ]]; then
  print Dynamic mode activated, batch size is $batch
  dynamic
  total_logfiles=$(ls *.dlg | wc -w)
  successful_jobs=$(grep "Successful\|Writing output ... done\|Total elapsed time" *.dlg | wc -l)
  print "Done; $total_logfiles docking logfiles started, $successful_jobs successfully completed." 
  
  print "ligID: $ligid  Date: $(date)"
  if (( successful_jobs < total_logfiles )); then
    print "$((total_logfiles-successful_jobs)) appear to have failed. Run rankAD.ksh -c to check for errors."
  fi
else
  print One-off mode activated
  oneoff
  sleep 11
fi
}


function distdock {
if [[ ! -e hostnames.txt ]]; then
  print "$hostnames" > hostnames.txt
fi
rm -f docklist_* autodist_*tmp
mode_selection=y
autodist -d
}


function vina {
print Vina function activated
if [[ ! -e hostnames.txt ]]; then
  print "$hostnames" > hostnames.txt
fi

if [[ ! -e vina.dpf ]]; then
  error="vina.dpf not found; please run autoAD.ksh -j to generate it."
  errhan
fi
mode_selection=z
rm -f docklist_* autodist_*tmp xscore.log
autodist -d
}


function distDOCK6 {
if [[ ! -e hostnames.txt ]]; then
  print "$hostnames" > hostnames.txt
fi
rm -f docklist_* autodist_*tmp
mode_selection=f
autodist -d
}


function run_xscore {
export XTOOL_HOME=/usr/people/douglas/programs/xscore_v1.2.1/
export XTOOL_PARAMETER=$XTOOL_HOME/parameter
export XSCORE_PARAMETER=$XTOOL_HOME/parameter
export XTOOL_BIN=$XTOOL_HOME/bin

xscore_protein=$1
xscore_pdbqt=$2
print "about to detect ${xscore_protein%.pdb?(qt)}_fixed.pdb" > xscore.log
if [[ ! -s ${xscore_protein%.pdb?(qt)}_fixed.pdb ]]; then
  print "${xscore_protein%.pdb?(qt)}_fixed.pdb does not exist; will create it from ${xscore_protein}." >> xscore.log
  /usr/people/douglas/programs/xscore_v1.2.1//bin/xscore -fixpdb ${xscore_protein} ${xscore_protein%.pdb?(qt)}_fixed.pdb >> xscore.log
fi 
print "finished detection" >> xscore.log

writetries=0
while [[ $(grep 'Writing output ... done.' ${xscore_pdbqt%.pdbqt}.dlg) != "Writing output ... done." ]]; do
  sleep 1
  ((writetries++))
  if ((writetries > 300)); then
    print "WARNING: $xscore_pdbqt on $(hostname -s) appears not to be finished (writetries are ${writetries}) - contents of Vinadlg is $(cat ${xscore_pdbqt%.pdbqt}.dlg)" >> xscore.log
    break
  fi
done
awk '/^MODEL 1$/,/^ENDMDL/' $xscore_pdbqt | grep "^[HA][ET][TO][AM][T ][M ]" | cut -c1-55 > model1_${xscore_pdbqt}.pdb; wait

writetries=0
while [[ ! -s model1_${xscore_pdbqt}.pdb ]]; do
  sleep 1
  ((writetries++))
  if ((writetries > 30)); then
    print "WARNING: awk did not extract $xscore_pdbqt to model1_${xscore_pdbqt}.pdb on $(hostname -s)" >> xscore.log
    print "Contents of model1_${xscore_pdbqt}.pdb is $(cat model1_${xscore_pdbqt}.pdb)" >> xscore.log
    print "Contents of $xscore_pdbqt is $(cat $xscore_pdbqt)" >> xscore.log
    break
  fi
done

/usr/local/bin/babel -ipdb model1_${xscore_pdbqt}.pdb -omol2 model1_${xscore_pdbqt}.mol2 -p 7.4 > /dev/null  2>&1 ; wait
writetries=0
while [[ ! -s model1_${xscore_pdbqt}.mol2 ]]; do
  sleep 1
  ((writetries++))
  if ((writetries > 30)); then
    print "WARNING: babel did not convert $xscore_pdbqt from model1_${xscore_pdbqt}.pdb to model1_${xscore_pdbqt}.mol2, contents of model1_${xscore_pdbqt}.pdb is $(cat model1_${xscore_pdbqt}.pdb)" >> xscore.log
    break
  fi
done


xscore_output=$(/usr/people/douglas/programs/xscore_v1.2.1//bin/xscore -fixmol2 model1_${xscore_pdbqt}.mol2 model1_${xscore_pdbqt}_fixed.mol2; wait)
writetries=0
while [[ ! -s model1_${xscore_pdbqt}_fixed.mol2 ]]; do
  if [[ $xscore_output != *This\ molecule\ is\ skipped* ]]; then
    sleep 1
    ((writetries++))
    if ((writetries > 30)); then
      cp model1_${xscore_pdbqt}.mol2 model1_${xscore_pdbqt}_fixed.mol2
      print "WARNING: fixmol skipped model1_${xscore_pdbqt}.mol2, writetries are ${writetries}, xscore output is $xscore_output" >> xscore.log
      break
    fi
  else
      cp model1_${xscore_pdbqt}.mol2 model1_${xscore_pdbqt}_fixed.mol2
      print "WARNING: fixmol skipped model1_${xscore_pdbqt}.mol2, writetries are ${writetries} (should be 0), xscore output is $xscore_output" >> xscore.log
      break
  fi
done

print "######################################################################
#                            XTOOL/SCORE                             # 
######################################################################
###
FUNCTION	SCORE
###
### set up input and output files ------------------------------------
###
#
RECEPTOR_PDB_FILE    ./${xscore_protein%.pdb?(qt)}_fixed.pdb
#REFERENCE_MOL2_FILE  none
#COFACTOR_MOL2_FILE  none 
LIGAND_MOL2_FILE     ./model1_${xscore_pdbqt}_fixed.mol2
#
OUTPUT_TABLE_FILE    ./xscore_${xscore_pdbqt}.table
OUTPUT_LOG_FILE      ./xscore_${xscore_pdbqt}.log
###
### how many top hits to extract from the LIGAND_MOL2_FILE?
###
NUMBER_OF_HITS       0 
HITS_DIRECTORY       ./hits.mdb 
###
### want to include atomic binding scores in the resulting Mol2 files?
###
SHOW_ATOM_BIND_SCORE	NO		[YES/NO]
###
### set up scoring functions -----------------------------------------
###
APPLY_HPSCORE         NO             	[YES/NO]
	HPSCORE_CVDW  0.004 
	HPSCORE_CHB   0.053
	HPSCORE_CHP   0.011
	HPSCORE_CRT  -0.061
	HPSCORE_C0    3.448
APPLY_HMSCORE         YES             	[YES/NO]
	HMSCORE_CVDW  0.004
	HMSCORE_CHB   0.094
	HMSCORE_CHM   0.394
	HMSCORE_CRT  -0.099
	HMSCORE_C0    3.585
APPLY_HSSCORE         NO 	  	[YES/NO]
	HSSCORE_CVDW  0.004
	HSSCORE_CHB   0.069
	HSSCORE_CHS   0.004
	HSSCORE_CRT  -0.092
	HSSCORE_C0    3.349
###
### set up chemical rules for pre-screening ligand molecules ---------
###
APPLY_CHEMICAL_RULES    NO            [YES/NO]	
	MAXIMAL_MOLECULAR_WEIGHT      600.0
	MINIMAL_MOLECULAR_WEIGHT      200.0
	MAXIMAL_LOGP                  6.00
	MINIMAL_LOGP                  1.00
	MAXIMAL_HB_ATOM               8 
	MINIMAL_HB_ATOM               2 
###
###
" > xscore_${xscore_pdbqt}.input

xscore_output=$(/usr/people/douglas/programs/xscore_v1.2.1//bin/xscore ./xscore_${xscore_pdbqt}.input; wait)
writetries=0
while [[ ! -e xscore_${xscore_pdbqt}.table ]]; do
  sleep 1
  ((writetries++))
  if ((writetries > 30)); then
    print "WARNING: X-score did not run on $xscore_pdbqt on $(hostname -s), writetries are $writetries, output is $xscore_output" >> xscore.log
    break
  fi
done

xscore_score=$(awk '{print $6}' xscore_${xscore_pdbqt}.table | tail -1 | tr -d "[[:space:]]"); wait

#adds the score to only the first model in the Vina output file
sed -i "1,/REMARK VINA RESULT/ {/REMARK VINA RESULT:/a\
REMARK XSCORE RESULT:      $xscore_score
}" $xscore_pdbqt; wait

rm -f model1_${xscore_pdbqt}.pdb  model1_${xscore_pdbqt}_fixed.mol2 xscore_${xscore_pdbqt}.table xscore_${xscore_pdbqt}.input xscore_${xscore_pdbqt}.log; wait

}

function run_drugscore {

drugscore_protein=$1
drugscore_pdbqt=$2
binding_site=$(head -5 vina.dpf | tail -1)


if [[ ! -s calc_pocket.in ]] || [[ ! -s ${drugscore_protein%.pdb?(qt)}_noH_poc.pdb ]]; then 
  /usr/local/bin/babel -ipdb $drugscore_protein -opdb ${drugscore_protein%.pdb?(qt)}_noH.pdb -d 2>/dev/null
  if [[ ${binding_site##*.} != mol2 ]]; then
    /usr/local/bin/babel -ipdb $binding_site -omol2 ${binding_site%.*}.mol2 2>/dev/null
  fi
  print ${drugscore_protein%.pdb?(qt)}_noH.pdb ${binding_site%.*}.mol2 > calc_pocket.in
  ~douglas/programs/drugscore/SCRIPTS/calc_pocket calc_pocket.in 9 y
fi

awk '/MODEL 1/,/ENDMDL/' $drugscore_pdbqt | egrep "^[AH][TE][OT][MA][ T][ M]|MODEL|ENDMDL" | cut -c1-55 | /usr/local/bin/babel -ipdb -omol2 -d  2>/dev/null | sed "s/\(@<TRIPOS>MOLECULE\)/\n\1/" > docked.mol2; wait
drugscore=$(~douglas/programs/drugscore/SCRIPTS/drugscore PAIRSURF ${drugscore_protein%.pdb?(qt)}_noH.pdb docked.mol2 | tail -1 | awk '{print $5}'); wait
#adds the score to only the first model in the Vina output file
sed -i "1,/REMARK VINA RESULT/ {/REMARK VINA RESULT:/a\
REMARK DRUGSCORE RESULT:      $drugscore
}" $drugscore_pdbqt; wait

}


$sdf2pdbs 
$pdbs2pdbqts 
$preprecep
$makegpf
$makemaps
$makedpfs 
$vinadpf
$prepDOCK6 
$makeDOCK6dpf
$autodock
$vina
$dynvina
$vina_oneoff
$autodist
$listvina


