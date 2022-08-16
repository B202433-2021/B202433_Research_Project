#!/usr/bin/ksh

source /localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/bin/mglenv.sh > /dev/null 2>&1
#setenv ACHOME /usr/people/ED/s1706179/antechamber-1.27/
# set path=($ACHOME/bin $path)

# performs the dock6 docking and preparation. Whether to use MOPAC or Gasteiger charges can be specified.

function prepare_receptor {
#Deletes alternate side chains on input protein, adds hydrogens, creates ${1%.pdb}_fixed_H.pdb
#USAGE:$0 protein.pdb
  ../bindsite.ksh ${dir}_protein.pdb ${dir}_ligand.pdb 15
  grep 'ATOM' ${dir}_protein_site.pdb | grep 'ALA\|CYS\|ASP\|GLU\|PHE\|GLY\|HIS\|ILE\|LYS\|LEU\|MET\|ASN\|PRO\|GLN\|ARG\|SER\|THR\|VAL\|TRP\|TYR' > ${dir}_protein_trim.pdb
  #only ions supported by Vina/AD4
  awk '$NF ~ /CA/ { print $0 }' ${dir}_protein_site.pdb > ${dir}_tmp2.pdb
  awk '$NF ~ /MG/ { print $0 }' ${dir}_protein_site.pdb >> ${dir}_tmp2.pdb
  awk '$NF ~ /MN/ { print $0 }' ${dir}_protein_site.pdb >> ${dir}_tmp2.pdb
  awk '$NF ~ /ZN/ { print $0 }' ${dir}_protein_site.pdb >> ${dir}_tmp2.pdb
  grep 'ATOM' ${dir}_protein_site.pdb | grep 'HEM\|FE' >> ${dir}_tmp2.pdb

  cat ${dir}_tmp2.pdb >> ${dir}_protein_trim.pdb
  print 'TER' >> ${dir}_protein_trim.pdb 

print "Deleting and adding Hs with obabel to make ${1%.pdb}_fixed_H.pdb."
obabel -ipdb ${dir}_protein_trim.pdb -opdb -O ${dir}_protein_trim_fixed_H.pdb -p 7.4 2>/dev/null
rm -f tmp* 
python ../hydrogens.py ${dir}_protein_trim_fixed_H.pdb
sleep 3
}


function split_pdbqt_to_mops {
#Splits each of the docked ligands in vina_results.pdbqt or ligand_largestC.pdbqt, cats it and the protein to a MOPAC .mop input file for a 1SCF calculation.
#USAGE: $0 protein.pdb dockingresults.pdbqt

if [[ $2 == *_out.pdbqt ]] then
  #no_models=$(grep "^MODEL " $2 | wc -l)
  no_models=1
  print "Vina results detected; no_models is ${no_models}."
  print "Creating $no_models MOPAC input files."
  for (( c=1; c <= $no_models; c++ )); do
    awk '/^MODEL '$c'/,/^ENDMDL/' $2 | grep "^[AH][TE][OT][MA][ T][ M]" | cut -c1-55 |  sed "s/\(^[AH][TE][OT][MA][ T][ M].\{11\}\).\{3\}\(.*\)/\1LIG\2/" | obabel -ipdb -opdb -O tmp1.pdb -d  2> /dev/null
    obabel -ipdb tmp1.pdb -opdb -O ${2%.pdbqt}_${c}.pdb -p 7.4 2>/dev/null
    #rm -f tmp1* 
    print "Created ${2%.pdbqt}_${c}.pdb - this should contain a ligand."
    cat ${1%.pdb}_fixed_H.pdb ${2%.pdbqt}_${c}.pdb | grep "^[AH][TE][OT][MA][ T][ M]" > ${1%.*}_${2%.*}_${c}.pdb
    print "Created ${1%.*}_${2%.*}_${c}.pdb - this should contain the protein and the ligand (please check)."
    #optimise ligand
    #label=$(grep -m 1 "LIG" ${1%.*}_${2%.*}_${c}.pdb | cut -c22-27| tr -d ' ')
    #brc=")"
    #
    #  Change MOPAC settings
    #optimise ligand
    print "MOZYME PDB EPS=78.4 CUTOFF=5 PM6-D3H4 GNORM=20\n\n\n$(cat ${1%.*}_${2%.*}_${c}.pdb)" > ${1%.*}_${2%.*}_${c}.mop
    #print "1SCF MOZYME PDB EPS=78.4 PM6-D3H4 CUTOFF=5\n\n\n$(cat ${1%.*}_${2%.*}_${c}.pdb)" > ${1%.*}_${2%.*}_${c}.mop
    print "Created ${1%.*}_${2%.*}_${c}.mop MOPAC input file - this should contain the protein and the ligand, with MOPAC keywords."
  done 
  c=1
fi
}


function mopac_complex {
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/localdisk/software/mopac2016
/localdisk/software/mopac2016/MOPAC2016.exe $1
}


#calculates spheres
function prepsphere {

#10 A seems to be best
dist_cutoff=10

#Remove hydrogens - this is what the tutorial says to do for dms
print "Running obabel -pdb ${arglist[${argno}]%.*}.pdb -opdb -O ${arglist[${argno}]%.*}_noH.pdb -d"
obabel -ipdb ${1%.*}.pdb -opdb -O t.pdb -d
obabel -ipdb t.pdb -opdb -O ${1%.*}_noH.pdb -d
sleep 2 ; wc -l  ${1%.*}_noH.pdb ; sleep 2

print "Running /localdisk/software/dms/bin/dms  ${arglist[${argno}]%.*}_noH.pdb -n -w1.4 -v -o ${arglist[${argno}]%.*}_noH.ms"
/localdisk/software/dms/bin/dms ${1%.*}_noH.pdb -n -w1.4 -v -o ${1%.*}_noH.ms
rm -f INSPH OUTSPH temp* ${1%.*}_noH.sph
print "${1%.*}_noH.ms
R
X
0.0
4.0
1.4
${1%.*}_noH.sph" > INSPH; wait

print "Running /localdisk/software/sphgen/sphgen_cpp_threads"
#~douglas/programs/dock6/bin/sphgen 
/localdisk/software/sphgen/sphgen_cpp_threads
print "ls immediately before sphere selector"
ls -alt
print "Running ~/Research_Project/dock6/bin/sphere_selector ${1%.*}_noH_site.sph $2 $dist_cutoff"
~/Research_Project/dock6/bin/sphere_selector ${1%.*}_noH.sph $2 $dist_cutoff
mv selected_spheres.sph ${1%.*}_selected_spheres.sph

print "Y
$dist_cutoff
${1%.*}_selected_spheres.sph
1
${1%.*}_noH_box.pdb" > showbox.in; wait

print "Running ~/Research_Project/dock6/bin/showbox < showbox.in"
~/Research_Project/dock6/bin/showbox < showbox.in

print "${1%.*}_selected_spheres.sph
1
N
sphgen_cluster.pdb" > showsphere.in

print "Running ~/Research_Project/dock6/bin/showsphere < showsphere.in"
~/Research_Project/dock6/bin/showsphere < showsphere.in
print "Spheres finished"
}

#converts ligand and receptor into mol2 and adds MOPAC charges
function receptor_ligand {
print "Calculating receptor_ligand charges: $3"
if [[ $3 == MOPAC ]]; then
#1-_protein_trim.mol2
obabel -ipdb ${1%.pdb}_fixed_H.pdb -omol2 -O ${2%.pdb}_${3}_protein.mol2 2>/dev/null
obabel -ipdb ${2%.pdb}.pdb -omol2 -O ${2%.pdb}_${3}.mol2 2>/dev/null
typeset -i num=$(grep "@<TRIPOS>BOND" -B1 ${2%.pdb}_${3}_protein.mol2 | head -n 1| cut -c 4-8)
mopac_charges=($(awk '/NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS/,/NET CHARGE ON RESIDUES/' ${1%.pdb}_${2%.pdb}.out | grep -v DIPOLE | tail --lines=+3 | cut -c55-63))
protein_charges=("${mopac_charges[@]:0:${num}}")
ligand_charges=("${mopac_charges[@]:${num}}")

#print ligand charges
print=0
n=0
lineno=0
rm -f tmp.mol2
 while IFS= read -r line; do
      ((lineno++))
#  print .$line.  
      if (( $print == 1 )); then
#     print ."${mopac_charges[$n]}".
        printf "${line:0:68} % 0.4f\n"   "${ligand_charges[$n]}" >> tmp.mol2
        ((n++))
      else
        print "$line" >> tmp.mol2
      fi

      if [[ $line == '@<TRIPOS>ATOM' ]]; then
        print=1
      fi

      if (( n >= ${#ligand_charges[*]} )); then
        print=0
        n=0
      fi
    done <  ${2%.pdb}_${3}.mol2
mv tmp.mol2 ${2%.pdb}_${3}.mol2
 
#print protein charges  
print=0
n=0
lineno=0
rm -f tmp.mol2
 while IFS= read -r line; do
      ((lineno++))
      if (( $print == 1 )); then
        printf "${line:0:68} % 0.4f\n"   "${protein_charges[$n]}" >> tmp.mol2
        ((n++))
      else
        print "$line" >> tmp.mol2
      fi

      if [[ $line == '@<TRIPOS>ATOM' ]]; then
        print=1
      fi

      if (( n >= ${#protein_charges[*]} )); then
        print=0
        n=0
      fi
    done <  ${2%.pdb}_${3}_protein.mol2
mv tmp.mol2 ${2%.pdb}_${3}_protein.mol2
rm -f tmp.mol2

print "@<TRIPOS>SUBSTRUCTURE
  1      GARBAGE 1       PERM    0       ****    ****    0       ROOT" >> ${2%.pdb}_${3}_protein.mol2
edits=($(head -3 ${2%.pdb}_${3}_protein.mol2 | tail -1))
sed -i "s/${edits[*]}/ ${edits[0]} ${edits[1]} 1 ${edits[3]}/" ${2%.pdb}_${3}_protein.mol2
sleep 3
if [[ $(head -5 ${2%.pdb}_${3}_protein.mol2 | tail -1) == GASTEIGER ]]; then
  rm -f tmp
  print "$(head -4 ${2%.pdb}_${3}_protein.mol2)\n$(tail --lines=+6 ${2%.pdb}_${3}_protein.mol2)" >> tmp
  mv tmp ${2%.pdb}_${3}_protein.mol2
fi
sleep 2


# if prepare gasteiger charges
elif [[ $3 == gast ]]; then
#add charges to metal ions and gasteiger charges to everything else; !!! Fe set as +2 !!!
  failed=$(/localdisk/home/s2173344/Research_Project/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -A none -U none -r ${1%.pdb}_fixed_H.pdb | grep 'Sorry')
  for lines in ${failed[*]};do
  typeset -R2 atom=${lines}
  if [ ${atom} == ZN ];then
    sed -i 's/0.000 Zn/2.000 Zn/g' ${1%.pdb}_fixed_H.pdbqt
  elif [ ${atom} == MG ];then
    sed -i 's/0.000 Mg/2.000 Mg/g' ${1%.pdb}_fixed_H.pdbqt
  elif [ ${atom} == CA ];then
    sed -i 's/0.000 Ca/2.000 Ca/g' ${1%.pdb}_fixed_H.pdbqt
  elif [ ${atom} == FE ];then
    sed -i 's/0.000 Fe/2.000 Fe/g' ${1%.pdb}_fixed_H.pdbqt
  elif [ ${atom} == MN ];then
    sed -i 's/0.000 Mn/2.000 Mn/g' ${1%.pdb}_fixed_H.pdbqt
  fi
  done
  sleep 3
  obabel -ipdbqt ${1%.pdb}_fixed_H.pdbqt -omol2 -O ${2%.pdb}_${3}_protein.mol2
  sleep 5
  print "@<TRIPOS>SUBSTRUCTURE
  1      GARBAGE 1       PERM    0       ****    ****    0       ROOT" >> ${2%.pdb}_${3}_protein.mol2
edits=($(head -3 ${2%.pdb}_${3}_protein.mol2 | tail -1))
sed -i "s/${edits[*]}/ ${edits[0]} ${edits[1]} 1 ${edits[3]}/" ${2%.pdb}_${3}_protein.mol2
sleep 3
if [[ $(head -5 ${2%.pdb}_${3}_protein.mol2 | tail -1) == GASTEIGER ]]; then
  rm -f tmp
  print "$(head -4 ${2%.pdb}_${3}_protein.mol2)\n$(tail --lines=+6 ${2%.pdb}_${3}_protein.mol2)" >> tmp
  mv tmp ${2%.pdb}_${3}_protein.mol2
fi
sleep 2
dir=$4
#prepare ligand gast
obabel -ipdb ${2%.pdb}.pdb -omol2 -O ${2%.pdb}_${3}.mol2 2>/dev/null



#if bcc charges
else
charges=$(sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/p' ../dock6_6/${2%.pdb}_DOCKed_scored.mol2 | cut -c 69-77) 
typeset -F sum=0
for j in ${charges[@]};do
  ((sum=sum+j))
done
typeset -i round
if ((sum<0));then
  ((round=sum-.5))
else
  ((round=sum+.5))
fi
rm -f *.AC *.ACO *ac *am1bcc* ${dir}.out ${dir}.arc
print "AM1 MMOK GEO-OK PRECISE CHARGE=${round} T=10000 DUMP=10000\n\n\n$(cat ${dir}_ligand.pdb)" > ${dir}.mop
ssh -T mscoc1 "cd $(pwd);  setenv MOPAC_LICENSE /opt/mopac; setenv LD_LIBRARY_PATH /opt/mopac;  /opt/mopac/MOPAC2016.exe ${dir}.mop "
#while [[ $(tail -1 ${dir}.out) != " == MOPAC DONE ==" ]]; do
#  sleep 10
#done
obabel -imopout ${dir}.out -omol2 > ${dir}_mopac.mol2
sleep 1

print=0
n=0
lineno=0
mopac_charges=($(awk '/NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS/,/NET CHARGE ON RESIDUES/' ${dir}.out | grep -v DIPOLE | tail --lines=+3 | cut -c55-63))
rm -f tmp.mol2
 while IFS= read -r line; do
      ((lineno++))
#  print .$line.  
      if (( $print == 1 )); then
#     print ."${mopac_charges[$n]}".
        printf "${line:0:68} % 0.4f\n"   "${mopac_charges[$n]}" >> tmp.mol2
        ((n++))
      else
        print "$line" >> tmp.mol2
      fi

      if [[ $line == '@<TRIPOS>ATOM' ]]; then
        print=1
      fi

      if (( n >= ${#mopac_charges[*]} )); then
        print=0
        n=0
      fi
    done <  ${dir}_mopac.mol2
mv tmp.mol2 ${dir}_mopac.mol2
ssh -T mscoc1 "cd $(pwd);  setenv ACHOME /usr/people/ED/s1706179/antechamber-1.27; set path=($ACHOME/bin $path); ~/antechamber-1.27/exe/antechamber -i ${dir}_mopac.mol2 -fi mol2 -fo ac -o ${dir}.ac; sleep 1; ~/antechamber-1.27/exe/am1bcc -i ${dir}.ac  -f ac -o ${dir}_am1bcc.ac -j 5"
sleep 5
am1_charges=($(grep 'ATOM' ${dir}_am1bcc.ac | cut -c 56-65))
while IFS= read -r line; do
      ((lineno++))
#  print .$line.  
      if (( $print == 1 )); then
#     print ."${am1_charges[$n]}".
        printf "${line:0:68} % 0.4f\n"   "${am1_charges[$n]}" >> tmp.mol2
        ((n++))
      else
        print "$line" >> tmp.mol2
      fi

      if [[ $line == '@<TRIPOS>ATOM' ]]; then
        print=1
      fi

      if (( n >= ${#am1_charges[*]} )); then
        print=0
        n=0
      fi
    done <  ${dir}_mopac.mol2
mv tmp.mol2 ${2%.pdb}.mol2


fi
}

#calculates grid
function grid {

print "compute_grids               yes
grid_spacing                   0.4
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             9
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  ${2%.pdb}_${3}_protein.mol2
box_file                       ${1%.pdb}_fixed_H_noH_box.pdb
vdw_definition_file            /localdisk/home/s2173344/Research_Project/dock6/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              ${2%.pdb}_${3}_grid" > ${2%.pdb}_${3}_grid.in; wait


print "Running ~/Research_Project/dock6/bin/grid -i ${2%.pdb}_${3}_grid.in"
~/Research_Project/dock6/bin/grid -i ${2%.pdb}_${3}_grid.in 2>/dev/null
}

function makeDOCK6dpf {

print "conformer_search_type                                        flex
write_fragment_libraries                                     no
ligand_atom_file                                             ${2%.pdb}_${3}.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ${1%.pdb}_fixed_H_selected_spheres.sph
max_orientations                                             1000
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
use_clash_overlap                                            no
clash_overlap                                                0.5
write_growth_tree                                            no
bump_filter                                                  no
bump_grid_prefix                                             ${2%.pdb}_${3}_grid
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
grid_score_grid_prefix                                       ${2%.pdb}_${3}_grid
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
continuous_score_primary                                     no
cont_score_rec_filename                                      ${2%.pdb}_${3}_protein.mol2
cont_score_dielectric                                        4
cont_score_att_exp                                           6
cont_score_rep_exp                                           12
cont_score_rep_rad_scale                                     1.0
cont_score_use_dist_dep_dielectric                           yes
cont_score_vdw_scale                                         1.0
cont_score_es_scale                                          1.0
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
vdw_defn_file                                                /localdisk/home/s2173344/Research_Project/dock6/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               /localdisk/home/s2173344/Research_Project/dock6/parameters/flex.defn
flex_drive_file                                              /localdisk/home/s2173344/Research_Project/dock6/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ${2%.pdb}_${3}_DOCKed
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
multigrid_score_secondary                                    no
SASA_score_secondary                                         no" > ${2%.pdb}_${3}.dockin
sleep 2
~/Research_Project/dock6/bin/dock6 -i ${2%.pdb}_${3}.dockin > ${2%.pdb}_${3}_DOCKed.dlg

# convert dock output files to pdb format (for rmsd calculation) and copy to appropriate directory
obabel -imol2 ${2%.pdb}_${3}_DOCKed_scored.mol2 -opdb -O ${2%.pdb}_${3}_DOCKed_scored.pdb
cp ${2%.pdb}_${3}_DOCKed_scored.pdb ../${3}/${2%.pdb}_${3}_DOCKed_scored.pdb
}

function run {
dir=$1
x=$2
charge=$3
#make grid
ls -alt
grid ${dir}_protein_trim.pdb ${dir}_ligand1_out_${x}.pdb ${charge}
sleep 2
#make dpf and dock
ls -alt
makeDOCK6dpf ${dir}_protein_trim.pdb ${dir}_ligand1_out_${x}.pdb ${charge}
print "makeDOCK6dpf just run"
ls -alt
}

function mopac {
  
      unset out
      for (( c=1; c <=${no_models}; c++ )); do
        mopac_complex ${1%.*}_${2%.*}_${c}.mop &
        sleep 5
        print $c away
        out+=($c)
    done 
    print all away


    #Check on the progress of the MOPAC jobs
    while (( ${#out[*]} > 0 )); do
      sleep 5
      for file in ${out[*]}; do
        if [[ $(tail -1 ${1%.*}_${2%.*}_${file}.out) == " == MOPAC DONE ==" ]]; then
          print "MOPAC finished on ${1%.*}_${2%.*}_${file}.out"
          print removing $file from outstanding
          ((uset=file-1))
          unset out[$uset]
          print "outstanding is now ${out[*]}"
        fi
    done
  done
  print "All MOPAC jobs finished."
  sleep 5
}

function dock_all {
dir=$1 
    cd $dir
    # get rid of existing stuff
    # rm *out* *trim* *.txt *.in sphgen_cluster.pdb tmp1.pdb t.pdb *_ligand.pdb *_protein_site.pdb  *_tmp2.pdb  INSPH
    obabel -imol2 ${dir}_ligand.mol2 -opdb -O ${dir}_ligand.pdb  # making sure correct input files (not sure if have to do if vina)
    # obabel -imol2 ${dir}_ligand.mol2 -opdbqt -O ${dir}_ligand1_out.pdbqt # making sure correct input files (not sure if have to do if vina)
    # cp ${dir}_ligand.mol2 ${dir}_ligand1_out_1.mol2 # making sure correct input files (not sure if have to do if vina)
    ls -alt
    # no_models=$(grep "^MODEL " ${dir}_ligand1_out.pdbqt | wc -l)
    no_models=1 # only do for top vina pose 
    # print "no_models = ${no_models}"
    prepare_receptor ${dir}
    ls -alt
    split_pdbqt_to_mops ${dir}_protein_trim.pdb ${dir}_ligand1_out.pdbqt
    ls -alt
    #calculate MOPAC charges
    # mopac ${dir}_protein_trim.pdb ${dir}_ligand1_out.pdbqt
    ls -alt 
    #calculate spheres
    obabel -ipdb ${dir}_ligand1_out_1.pdb -omol2 -O ${dir}_ligand1_out_1.mol2
    prepsphere ${dir}_protein_trim_fixed_H.pdb ${dir}_ligand1_out_1.mol2
    ls -alt
    rm -f *.dlg *.pdbqt
    ls -alt
    unset outstanding
    set -A charges gast
    for ((x=1;x<=${no_models};x++)); do
      for charge in ${charges[@]}; do
        print "Doing ${charge} receptor ligand charges"
        receptor_ligand ${dir}_protein_trim.pdb ${dir}_ligand1_out_${x}.pdb ${charge} ${dir}
        ls -alt
        sleep 2
        #dock: this is actually just a comment
        run ${dir} ${x} ${charge} &
        sleep 1
#        outstanding+=($x)
#        print "outstanding is now ${outstanding[*]}"
      done
    done

    sleep 10
    #check until finished
#    while (( ${#outstanding[*]} > 0 )); do 
#    sleep 30
 #     for file in ${outstanding[*]}; do
  #      if [[ $(tail ${dir}_ligand1_out_${file}_${charge}_DOCKed.dlg | grep "Total elapsed time") == *(?)+(Total elapsed time)*(?) ]]; then
   #       print removing $file from outstanding
    #      ((uset=file-1))
     #     unset outstanding[$uset]
     #     print "outstanding is now ${outstanding[*]}"
   #     fi
   #   done
  #  done   
    cd ..
    cp -R ${dir} MOPAC_docked_directories_backup/    
}
#main 
for dir in $(<$1); do
print $dir; sleep 1
    dock_all $dir &
done

