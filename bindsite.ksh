#!/usr/bin/ksh
###########MEASURES DISTANCE BETWEEN ATOMS IN 2 PDB FILES, WRITES OUT THOSE WITHIN X A
# USAGE="usage : $0 <protein.pdb> <ligand.pdb> <distcutoff in A>"
while [[ $1 == -* ]]
do
   case $1 in
      -a ) adtmacro=1 ;;
      -r ) residues=1 ;;
      *  ) error="This program takes no argument (got it?)." 
        errhan $0;;
   esac
   shift
done

# 1st 2 arguments must be readable regular files
if [[ (-f $1) && (-f $2) && (-r $1) && (-r $2) ]]
then :
else
    print "$USAGE"
    exit 2
fi

rm -f ${1%.*}_site.pdb
#grep -v  "^.\{16\}B" $1 > "$1"_tmp.$$
#cp  $1 "$1"_tmp.$$
egrep "^ATOM|^HETATM" $1 > ${1}_tmp.$$

n=0
o=0
#print o is $o and LC1 is $LC1 - ligand atom
#print n is $n and LC2 is $LC2 DIST is $DIST - protein atom
print "Measuring $3 A around ligand:"

COA1=($(egrep "^ATOM|^HETATM" $2 | cut -c31-38))
COA2=($(egrep "^ATOM|^HETATM" $2 | cut -c39-46))
COA3=($(egrep "^ATOM|^HETATM" $2 | cut -c47-54))

COB1=($(egrep "^ATOM|^HETATM" ${1}_tmp.$$ | cut -c31-38))
COB2=($(egrep "^ATOM|^HETATM" ${1}_tmp.$$ | cut -c39-46))
COB3=($(egrep "^ATOM|^HETATM" ${1}_tmp.$$ | cut -c47-54))

while (( $o < ${#COA1[*]} ))
do
  print "working on ligand atom $((o+1)) ${COA1[$o]}"
  while (( $n < ${#COB1[*]} ))
  do
    ((LE1=${COA1[$o]} - ${COB1[$n]}))
    ((LE2=${COA2[$o]} - ${COB2[$n]}))
    ((LE3=${COA3[$o]} - ${COB3[$n]}))

    ((DIST2 = LE1*LE1 + LE2*LE2 + LE3*LE3))
    DIST=$((sqrt($DIST2)))

    if (( $DIST <= $3 )); then
      writeline=$n; ((writeline++))
#      print writing out line $(head -$writeline "$1"_tmp.$$ | tail -1 ), o $o ${#COA1[*]} n $n ${#COB1[*]} DIST $DIST COA1 ${COA1[$o]} COB1 ${COB1[$n]} COA2 ${COA2[$o]} COB2 ${COB2[$n]}
       head -$writeline "$1"_tmp.$$ | tail -1 >> ${1%.*}_site.pdb
    fi

    ((n++))
  done
  n=0
  ((o++))
done
#done


siteatomlist=$(cut -c18-26 ${1%.*}_site.pdb | tr " " ":" | sort -u --key=1.6)
for siteatom in $siteatomlist
do
#  print siteatom ${siteatom//:/ }e
  if [[ -n $adtmacro ]]; then
    siteatom2=${siteatom:0:3}${siteatom:5:4}
    print "self.selectFromString(xor=False, log=False, res='${siteatom2// }', mols='', atoms='', intersect=False, negate=False, chains='', silent=True)" >> ${1%.pdb}_site.py_tmp.$$
  fi
  grep "^[HA][ET][TO][AM][T ][M ].\{11\}${siteatom//:/ }" "$1"_tmp.$$  | sort -u >> ${1%.pdb}_site.pdb_tmp.$$
done

mv ${1%.pdb}_site.pdb_tmp.$$ ${1%.pdb}_site.pdb
print Binding site has been written to ${1%.pdb}_site.pdb

print The following residues have been found to have at least 1 atom within $3 A of the ligand:
cut -c18-26 ${1%.pdb}_site.pdb | sort -u | sort -n --key=1.6


if [[ -n $adtmacro ]]; then
  mv ${1%.pdb}_site.py_tmp.$$ ${1%.pdb}_site.py
  print ADT macro has been written to ${1%.pdb}_site.py
fi
rm *.$$
