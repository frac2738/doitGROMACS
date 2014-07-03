#!/bin/bash
set -e   
###############################################################################
#     Version:      V1.0.3                                                    #
#     Last update:  02.07.14                                                  #
#     Author:       Francesco Carbone                                         #
#     Description:  Script to execute a bunch of stuff with gromacs           #
#     Updates :     - fix options                                             #
#                   - improved the help message                               #
#                   - improved the comments for all the functions             #
#                                                                             #
###############################################################################

################# Function declaration ################# 

helpMessage() {
   cat <<EOF

                     doitGROMACS.sh -  version 1.0.0  

Copyright (c) 2013-2014, University College London (UCL), Francesco Carbone

This script is designed to automatise the first step of a molecular dynamics
experiment (solvation and equilibration) and run some basic analyses on the
trajectory files (.xtc not .trr).
This script was written using GROMACS 4.6 and although it should work with any
previous versions, it is advise to check the commands before using a different
version.

Option   Type     Value       Description                  [PROBABLY REQUIRED]
--------------------------------------------------------------------------------
-[no]h   bool     yes         Print help info
-b       string   acrm        Set the location of gromacs binaries
                              acrm     -> Darwin building computer
                              emerald  -> Emerald cluster
                              bear     -> Personal laptop
-n       int      wt          Set the name
-t       int      200         Set the simulation length
-r       string   r1          Set the number of the replica
-k       int      400         Set the temperature in KELVIN


Option   Type     Value       Description      [OPTIONAL (function dependant)]
--------------------------------------------------------------------------------
-s       string   .tpr        .tpr file          
-f       string   .xtc        trajectory file
-c       string   .pdb        pdb file to use to start a simulation
-e       string   .edr        Energy file 

NOTE: In my simualtions all the output are printed in this format:
                           NAME_rX_TIME
      where NAME is the name of the mutation (306r), rX is the replica number
     (r1,r2,...) and TIME is the simulation time. As a consequence this script 
      takes and process outputs names in this form.
         

EOF
}


# modVim: it takes a text file and it converts all the comment character to "#" using vim (or vi).
# This is to avoid the use of "@" as formatting character in grace files.
modVim() {
   ex $1 << EOEX
      :%s/@/#/g
      :wq
EOEX
}

# inputs: it takes a pdb file and generate all the outputs required for the energy
# minimisation; it acts in two steps:
# 1) it creates a topology file depending on the force field (ff) and the water model chosen;
# 2) it creates a box around the protein and it solvates it.
# REQUIRED FILES:
# OUTPUTS:

inputs() { 
   # create the topology.
   # 6 = AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
   # 1 = TIP3P
   (echo "6"; echo "1") | $path/pdb2gmx -f $pdb1 -o $name1"_processed.pdb" -p topology.top -ignh -v

   read -e -p "select the type of box? " boxtype
   read -e -p "select the protein-box distance? (nm) " distedge
   $path/editconf -f $name1"_processed.pdb" -o $name1"_inbox.pdb" -bt $boxtype \
      -d $distedge -c
   $path/genbox -cp $name1"_inbox.pdb" -cs spc216.gro -o $name1"_sol.pdb"      \
      -p topology.top
   $path/grompp -f G6PD_ionization.mdp -c $name1"_sol.pdb" -p topology.top     \
      -o input_ioni.tpr 
   read -e -p "how many ions do you need? " pioni
   $path/genion -s input_ioni.tpr -p topology.top -o $name1"_ioni.pdb"         \
      -pname NA -nname CL -np $pioni
   # for negative ions use the -np flag
}

# energy_minimization: function that minimises a structure.
# REQUIRED FILES:
# OUTPUTS: 
energy_minimization() {
   $path/grompp -f $temp1"_min.mdp" -c $name1"_ioni.pdb" -p topology.top       \
      -o input_min.tpr
   $path/mdrun -s input_min.tpr -deffnm $name1"_min" -v  
   # export the potential energy profile
   echo Potential | $path/g_energy -f $name1"_min.edr" -o $name1"_potential.xvg"
   modVim $name1"_potential.xvg"
   # plot the potential energy profile using ggplot
   #GGplot $name1"_potential.xvg"
}  

nvt() {
   $path/grompp -f $temp1"_nvt.mdp" -c $name1"_min.gro" -p topology.top        \
      -o input_nvt.tpr
   $path/mdrun -s input_nvt.tpr -deffnm $name1"_nvt" -v
   # export the temperature profile 
   echo Temperature | $path/g_energy -f $name1"_nvt.edr" -o $name1"_temperature.xvg"
   modVim $name1"_temperature.xvg"
   # plot the temperature profile using ggplot
} 

npt() {
   $path/grompp -f $temp1"_npt.mdp" -c $name1"_nvt.gro" -p topology.top        \
      -o input_npt.tpr
   $path/mdrun -s input_npt.tpr -deffnm $name1"_npt" -v  
   # export the pressure profile 
   echo Pressure | $path/g_energy -f $name1"_npt.edr" -o $name1"_pressure.xvg"
   modVim $name1"_pressure.xvg"
   # export the density profile
   echo Density | $path/g_energy -f $name1"_npt.edr" -o $name1"_density.xvg"
   modVim $name1"_density.xvg"
} 

sim_conditions() {
   # export potential energy, temperature, pressure and density profiles 
   echo Potential | $path/g_energy -f $energy -o $nameprod"_potential.xvg"
   modVim $nameprod"_potential.xvg"
   echo Temperature | $path/g_energy -f $energy -o $nameprod"_temperature.xvg"
   modVim $nameprod"_temperature.xvg"
   echo Pressure | $path/g_energy -f $energy -o $nameprod"_pressure.xvg"
   modVim $nameprod"_pressure.xvg"
   echo Density | $path/g_energy -f $energy -o $nameprod"_density.xvg"
   modVim $nameprod"_density.xvg"
   # plot using ggplot
}

clean_trj() {
   # removing water molecules from the trajectory file
   echo Protein | $path/trjconv -s "input_"$timens".tpr" -f $nameprod".xtc"    \
      -o $nameprod"_only.xtc"
   # removing water molecules from the tpr file 
   echo Protein | $path/tpbconv -s "input_"$timens".tpr" -o $nameprod"_only.tpr"
   # creating a clean .gro file 
   echo Protein | $path/trjconv -s $nameprod"_only.tpr" -f $nameprod"_only.xtc"\
      -o $nameprod"_only.gro" -dump 0
   # account for the periodicity (nojump)
   echo Protein | $path/trjconv -s $nameprod"_only.tpr" -f $nameprod"_only.xtc"\
      -o $nameprod"_nojump.xtc" -pbc nojump
   # account for the periodicity (fitting) 
   (echo "Backbone"; echo "Protein") | $path/trjconv -s $nameprod"_only.tpr"   \
      -f $nameprod"_nojump.xtc" -o $nameprod"_fit.xtc" -fit rot+trans
   # remove intermediate files and rename them in less complicated way
   rm $nameprod"_nojump.xtc" $nameprod".xtc" $nameprod"_only.xtc"
   mv $nameprod"_fit.xtc" $nameprod".xtc"
   mv $nameprod"_only.tpr" $nameprod".tpr"
   mv $nameprod"_only.gro" $nameprod".gro"
}

rmsdf() {
   # calculating the RMSD 
   (echo "Backbone"; echo "Backbone") | $path/g_rms -s $nameprod".tpr"         \
      -f $nameprod".xtc" -o $nameprod"_rmsd.xvg" -tu ns
   # calculating the radius of Gyration 
   echo Protein | $path/g_gyrate -s $nameprod".tpr" -f $nameprod".xtc"         \
      -o $nameprod"_rgyration.xvg"
   echo " calculate rmsf for the backbone"
   echo Backbone | $path/g_rmsf -s $nameprod".tpr" -f $nameprod".xtc"          \
      -o $nameprod"_rmsf_bb.xvg" -oq $nameprod"_rmsf_bb.pdb" -res 
      # res=averages for each residues
   echo "calculate rmsf for the sidechains"
   echo SideChain | $path/g_rmsf -s $nameprod".tpr" -f $nameprod".xtc"         \
      -o $nameprod"_rmsf_sc.xvg" -oq $nameprod"_rmsf_sc.pdb" -res 
}

cluster_analysis() {
   
   read -e -p " skip? " skip
   if [ ! -d ./clusters_$nameprod2 ] ; then
      mkdir clusters_$nameprod2
   fi
   cd clusters_$nameprod2
   # check if the rmsd matrix exists
   while [ ! -f ./rmsd-matrix.xpm ]
   do
      # create the rmsd-matrix (using the backbone for the calculation)
      (echo "Backbone"; echo "Backbone") | $path/g_rms -s ../$tpr -f ../$trj   \
         -f2 ../$trj -m rmsd-matrix.xpm -dist rmsd-distribution.xvg -tu ns -skip $skip
      # improve rmsd matrix plot
      $path/xpm2ps -f rmsd-matrix.xpm -o rmsd-matrix.eps
   done
   #check the distribution file and decide the cutoff
   xmgrace rmsd-distribution.xvg
   read -e -p "Which cutoff do you want to use? " cutoff
   method='gromos'
   # cluster analysis on the rmsd matrix (Using the backbone for the calculation)
   (echo "Backbone"; echo "Backbone") | $path/g_cluster -s ../$tpr -f ../$trj  \
      -dm rmsd-matrix.xpm -o clusters.xpm -sz clusters-size.xvg                \
      -clid clusters-ovt.xvg -cl clusters.pdb -cutoff $cutoff -method $method  \
      -tu ns -skip $skip
   # to visualize in pymol use
   # split_states clusters
   # delete clusters
   # dss
   # show cartoon
   cd ..
}

pca() {
   read -e -p " dt? " dt
   if [ ! -d ./PCA_$nameprod2 ] ; then
      mkdir PCA_$nameprod2
   fi
   cd PCA_$nameprod2
   # check if the covariance matrix exists
   while [ ! -f ./covariance.xpm ]
   do
      # calculating the covariance matrix (C-alpha), eigenvalues and eigenvectors
      (echo "C-alpha"; echo "C-alpha") | $path/g_covar -s ../$tpr -f ../$trj   \
         -o eigenvalues.xvg -v eigenvectors.trr -xpma covariance.xpm -tu ns -dt $dt
      # improve covariance matrix ploting
      $path/xpm2ps -f covariance.xpm -o covariance.eps
   done
   xmgrace eigenvalues.xvg
   # determinare su quali autovalori lavorare
   read -e -p "how many eigenvalues do you want to use for the analysis? " range
   (echo "C-alpha"; echo "C-alpha") | $path/g_anaeig -v eigenvectors.trr       \
      -s ../$tpr -f ../$trj -first 1 -last $range -proj "projection-1"$range".xvg" -tu ns
   for ((i = 1; i <= range; i=i+1))
   do
      (echo "C-alpha"; echo "C-alpha") | $path/g_anaeig -v eigenvectors.trr    \
         -s ../$tpr -f ../$trj -first $i -last $i -nframes 100                 \
         -extr "ev"$i".pdb" -filt "ev"$i".xtc" -proj "ev"$i".xvg" -tu ns
      $path/g_analyze -f "ev"$i".xvg" -cc "ev"$i"-cc.xvg" -ac "ev"$i"-ac.xvg"
   done
   # calculate 2d projections and FES
   (echo "C-alpha"; echo "C-alpha") | $path/g_anaeig -v eigenvectors.trr       \
      -s ../$tpr -f ../$trj -first 1 -last 2 -2d 2d-12.xvg
   $path/g_sham -f 2d-12.xvg -ls gibbs-12.xpm -notime
   $path/xpm2ps -f gibbs-12.xpm -o gibbs-12.eps -rainbow red
   cd ..
}

GGplot() {
   write something
}

# deprecated will be rewritten fro SAS analysis
patches() {
   # create the pdb directory
   if [ ! -d ./patches_$nameprod2 ] ; then
      mkdir patches_$nameprod2
   fi
   cd patches_$nameprod2
   read -e -p "how long is your simulation? (ps) " simps
   read -e -p "jump? (ps)" jump
   for (( i = 0; i<= simps ; i=i+jump ))
   do
      echo Protein | $path/trjconv -s ../$tpr -f ../$trj -o $i".pdb" -dump $i
   done
   # run Tom's patches / write perl and call it here!!!
   /home/bsm3/zcbtfo4/public/automatic_patches -r 8 -t normal -l \
      $nameprod"_patchlog" -o $nameprod"_patches" -
   cd ..
}

################# end function declaration ################# 

################# The program begins here ################

while getopts ":hb:n:t:r:k:s:f:c:" opt; do
   case $opt in
      h) helpMessage; exit    ;;
      b) cpu=$OPTARG          ;;
      n) name1=$OPTARG        ;;
      t) timens=$OPTARG       ;;
      r) replica1=$OPTARG     ;;
      k) temp=$OPTARG         ;;
      s) tpr=$OPTARG          ;;
      f) trj=$OPTARG          ;;
      c) pdb1=$OPTARG         ;;   
      e) energy=$OPTARG       ;;
      \?) helpMessage;  exit  ;;
   esac
done

# check if no arguments are passed and in that case print the help message.
if ( ! getopts ":hb:n:t:r:k:s:f:c:e:" opt); then
	helpMessage;   exit
fi

# depending on the machine the script is running, locate both gromacs and R executables.
case $cpu in
   acrm) path='/acrm/usr/local/apps/gromacs/bin'
         Rpath='/export/francesco/R-3.1.0/bin/R'  ;;
   emerald) path='/apps/gromacs/4.6.3/bin'   ;;
   bear) path='/usr/local/gromacs/bin'
         Rpath='/usr/local/bin/R'  ;;
   *) echo " ERROR!! ERROR!! no GROMACS executable found  ERROR!! ERROR!! "
esac 

########################

echo " ----------------- THIS ARE YOUR OPTIONS ----------------- "
echo " 1 -  Starting from scratch "
echo " 2 -  Starting from E-minimization "
echo " 3 -  Starting from NVT "
echo " 4 -  Starting from NPT "
echo " 5 -  Remove water from a trajectory file "
echo " 6 -  Plot all the simulation conditions (U-T-P-density)"
echo " 7 -  Calculate RMSD, GYRATION RADIUS and RMSF [for backbone and sidechains] "
echo " 8 -  Cluster analysis "
echo " 9 -  PCA analysis "
echo " 10 - Do 8 and 9 "
echo " 11 - Patch analysis [soon] "
echo " 12 - Plot [soon] "
echo "-----------------------------------------------------------"


read -e -p "What do you want to do? " choice
case $choice in
   1|2|3|4|5|6|7|8|9|10|11|12)
   ;; 
   *)
   echo " ERROR!! ERROR!! $choice not listed ERROR!! ERROR!! "
   exit
   ;;
esac

# Check optional dependencies [arguments required]
case $choice in
   1) if [ -z "${pdb1+x}" ]; then
      ls
      read -e -p "which pdb do you want to use ? " pdb1
   fi;;
   6) if [ -z "${energy+x}" ]; then
      read -e -p "which edr do you want to use? " energy
      fi;;
   5|7|8|9|10)
      if [ -z "${tpr+x}" ] || [ -z "${trj+x}" ]; then
         ls
         read -e -p "which tpr do you want to use? " tpr
         read -e -p "which trj do you want to use? " trj
      fi;;
esac

nameprod="${name1}_${replica1}_${timens}"
nameprod2="${name1}_${replica1}"

case $choice in
   1 )
      inputs && energy_minimization && nvt && npt  ;;
   2)
      energy_minimization && nvt && npt   ;;
   3)
      nvt && npt  ;;
   4)
      npt   ;;
   5)
      clean_trj   ;;
   6) 
      sim_conditions ;;
   7)
      rmsdf ;;
   8)
      i=1
      cluster_analysis
      read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
      while [ "$ramen" == "yes" ] 
      do
         mv clusters_$nameprod2 clusters"$i"_$nameprod2
         cluster_analysis
         i++
         read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
      done  ;;
   9)
      pca   ;;
   10)
      cluster_analysis && pca ;;
   11)
      patches  ;;
   12)
      plot  ;;
esac
 
################ The program ends here ################
