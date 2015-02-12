#!/bin/bash
set -e   
###############################################################################
#     Version:      V1.0.4                                                    #
#     Last update:  17.10.14                                                  #
#     Author:       Francesco Carbone, UCL                                    #
#     Description:  Script to execute a bunch of stuff with gromacs           #
#     Usage:        Execute   "doitGROMACS_wopts.sh -h"                       #
#     Updates:      - fix a bug in the edr function                           #
#                   - add sas function                                        #
###############################################################################

#############
helpMessage() {
   cat <<EOF

                        doitGROMACS.sh -  v 1.0.4  

   Copyright (c) 2013-2014, University College London (UCL), Francesco Carbone

   This script is designed to automatise the first step of a molecular dynamics
   experiment ( solvation and equilibration) and run some basic analyses on the
   trajectory files (.xtc), using GROMACS tools.
   This script  was written using GROMACS 4.6  and although it should work with
   any previous versions, it is advise  to check the  validity of the  commands
   before  using  a different version.


   Option   Type     Value       Description                  [ALWAYS REQUIRED]
   ----------------------------------------------------------------------------
   -[no]h   bool     yes         Print help info
   -b       string   acrm        Set the location of binaries
                                 acrm     -> Darwin building computer
                                 emerald  -> Emerald cluster
                                 bear     -> Personal laptop
   -n       int      wt          Set the name


   Option   Type     Value       Description    [OPTIONAL (function dependant)]
   ----------------------------------------------------------------------------
   -t       int      200         Set the simulation length
   -k       int      400         Set the temperature in KELVIN
   -s       string   .tpr        .tpr file          
   -f       string   .xtc        Trajectory file
   -c       string   .pdb        Pdb file to use to start a simulation
   -e       string   .edr        Energy file 
   -u       string   txt file    Unres simulation file


   NOTE 1:  The "b" flag is used to set both gromacs and R binaries, so this 
            script MUST BE EDITED depending on the machine used.

   NOTE 2:  In my simualtions all the output are printed in this format:
                              NAMErX_TIME
            where  NAME is the name of the mutation ( 306r ), rX is the replica
            number (r1,r2,...) and TIME is the simulation time.As a consequence
            this script takes and process outputs names in this form.
            
EOF
}

#############
doitOptions() {
   cat <<EOF

                             -----------
   ------------------------- doitOPTIONS -------------------------
                             -----------

         1  - Starting from scratch ** 
         2  - Starting from E-minimization ** 
         3  - Starting from NVT **
         4  - Starting from NPT **
         5  - Remove water from a trajectory file 
         6  - Plot all the simulation conditions (U-T-P-density)
         7  - Calculate RMSD, GYRATION RADIUS and RMSF [backbone & sidechains] 
         8  - Cluster analysis 
         9  - PCA analysis 
         10 - Cluster and PCA
         11 - SAS analysis  
         12 - Plot with ggplot(R) ## [not yet implemented] 

   ----------------------------------------------------------------------------

   ** Options that require a parameter file (.mdp), NOT included and that MUST 
      be placed in the same directory; the functions only accept these files:
                        
                        TEMP_min.mdp   
                        TEMP_nvt.mdp
                        TEMP_npt.mdp
                        TEMP_md.mdp

      with TEMP = temperature in Kelvin (310, 400, ...)

   ## Option that requires a R scripts. The script is NOT included and MUST be 
      passed through the variable "RscriptPATH" (see -h NOTE 1)
   

EOF
}
   
########
# DESCRIPTION: Replace all the "@" in a text file with "#" using sed. Usefull 
#              to avoid the use of "@" as formatting character in grace files.
modVim() {
   sed -i 's/@/#/g' "$1" # -i modify and save 
}

########
# DESCRIPTION:  It prepares all the input files required for a simulation. 
#               Using a force field and a a water model, it first generates a 
#               topology file (pdb2gmx), it creates a box around the protein 
#               (editconf), solfatate (genbox) and ionised (grompp) the system.
# REQUIREMENTS: .pdb + .mdp (gromacs parameter file)
# OUTPUTS:      
inputs() { 
   # if -c is not set prompt and exit
   if [ -z "${pdb1}" ]; then
      cat <<EOF
      
               #########################################
               #        a pdb is required (-c)         #
               #########################################
EOF
      helpMessage; exit ;
   fi
   
   # create the topology
   # 6 = AMBER99SB-ILDN protein, nucleic AMBER94 
   # 1 = TIP3P
   (echo "6"; echo "1") | $path/pdb2gmx -f $pdb1 -o $name1"_processed.pdb"     \
      -p topology.top -ignh -v
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
   # for negative ions use the -nn flag
}

#####################
# DESCRIPTION: It minimised a 
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

#####
# DESCRIPTION:
nvt() {
   $path/grompp -f $temp1"_nvt.mdp" -c $name1"_min.gro" -p topology.top        \
      -o input_nvt.tpr
   $path/mdrun -s input_nvt.tpr -deffnm $name1"_nvt" -v
   # export the temperature profile 
   echo Temperature | $path/g_energy -f $name1"_nvt.edr" -o $name1"_temperature.xvg"
   modVim $name1"_temperature.xvg"
   # plot the temperature profile using ggplot
} 

#####
# DESCRIPTION:
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

################
# DESCRIPTION:
sim_conditions() {
   # if -e is not set prompt and exit
   if [ -z "${energy}" ]; then
      cat <<EOF
      
               #########################################
               #    an energy file is required (-e)    #
               #########################################
EOF
      helpMessage; exit ;
   fi

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

###########
# DESCRIPTION:
clean_trj() {
   # if -s and -f are not set prompt and exit
   if [[ -z "${tpr+x}" && -z "${trj+x}" ]]; then
      cat <<EOF
      
               #########################################
               #     tpr and trj files required        #
               #########################################
EOF
      helpMessage; exit ;
   fi
   # removing water molecules from the trajectory file
   echo Protein | $path/trjconv -s $tpr -f $trj    \
      -o $nameprod"_only.xtc"
   # removing water molecules from the tpr file 
   echo Protein | $path/tpbconv -s $tpr -o $nameprod"_only.tpr"
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

#######
# DESCRIPTION:
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

#################
# DESCRIPTION:
clusterAnalysis() {
   # if -s and -f are not set prompt and exit
   if [[ -z "${tpr+x}" && -z "${trj+x}" ]]; then
      cat <<EOF
      
               #########################################
               #     tpr and trj files required        #
               #########################################
EOF
      helpMessage; exit ;
   fi
   
   read -e -p " skip? " skip
   if [ ! -d ./clusters_$name1 ] ; then
      mkdir clusters_$name1
   fi
   cd clusters_$name1
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

repeatClusterAnalysis() {
   mv clusters_$name1 clusters"$i"_$name1
   clusterAnalysis
   i++
   read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
}

#####
# DESCRIPTION:
pca() {
   # if -s and -f are not set prompt and exit
   if [[ -z "${tpr+x}" && -z "${trj+x}" ]]; then
      cat <<EOF
      
               #########################################
               #     tpr and trj files required        #
               #########################################
EOF
      helpMessage; exit ;
   fi
   
   read -e -p " dt? " dt
   if [ ! -d ./PCA_$name1 ] ; then
      mkdir PCA_$name1
   fi
   cd PCA_$name1
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

#############
# DESCRIPTION:
sasAnalysis() {
   # if -s and -f are not set prompt and exit
   if [[ -z "${tpr+x}" && -z "${trj+x}" ]]; then
      cat <<EOF
      
               #########################################
               #     tpr and trj files required        #
               #########################################
EOF
      helpMessage; exit ;
   fi
   
   # create the directory
   if [ ! -d ./sas_$name1 ] ; then
      mkdir sas_$name1
   fi
   cd sas_$name1
   (echo "Protein"; echo "Protein") | $path/g_sas -s ../$tpr -f ../$trj        \
      -o $name1"_area.xvg" -tv $name1"_volume.xvg" -q $name1"_connelly.pdb"    \
      -dt 500 -probe 1.6 
   # probe 1.6 nm perchè 16A
   cd ..
}

########
# DESCRIPTION:
GGplot() {
   # check if the files esistono e nel caso run simconditions x crearli
   write something
}

############
# DESCRIPTION: if flag -u exists take the file and run unres analyses.pl
#              controllare all'inizio l'esistenza di -u e in quel caso -b -n ..
               # non sono necessari
unresAnalyses() {
   echo $unres
}

###############################################################################

############################ The program begins here ##########################

while getopts "hu:b:n:t:k:s:f:c:e:" opt; do
   case $opt in
      h) helpMessage; exit    ;;
      u) unres=$OPTARG        ;;
      b) cpu=$OPTARG          ;;
      n) name1=$OPTARG        ;;
      t) timens=$OPTARG       ;;
      k) temp1=$OPTARG        ;;
      s) tpr=$OPTARG          ;;
      f) trj=$OPTARG          ;;
      c) pdb1=$OPTARG         ;;   
      e) energy=$OPTARG       ;;
      \?) helpMessage;  exit  ;;
   esac
done

checkFlags() {
   # check -b flag 
   if [[ -z "${cpu}" && -z "${name1}" ]]; then
      cat <<EOF

               #########################################
               #  Both flags -b AND -n are required    #
               #########################################
EOF
      helpMessage; exit;
   fi 
}



# check se -u è presente e in quel caso chiamare unres, altrimenti eseguire 
# gromacs

if [ -n "${unres}" ]; then
   unresAnalyses; exit;   
else
   #if ( ! getopts ":hb:n:t:r:k:s:f:c:e:" opt); then
   #	helpMessage;   exit
   #fi

   checkFlags
   
   # Set the paths to the executables + location of all the R scripts.
   case $cpu in
      acrm) path='/acrm/usr/local/apps/gromacs/bin'
            Rpath='/export/francesco/R-3.1.0/bin/R'
            RscriptPATH='/export/francesco/Dropbox/scripts/Rscripts/Rplot_gromacs' ;;
      emerald) path='/apps/gromacs/4.6.3/bin'   ;;
      bear) path='/usr/local/gromacs/bin'
            Rpath='/usr/bin/R'  
            RscriptPATH='~/Dropbox/scripts/Rscripts/Rplot_gromacs';;
      *) echo "
                     ---- no executables found ---- 
         "
      helpMessage;   exit  ;;
   esac 

   # list the options 
   doitOptions

   # check the existance of the selected option and stop the execution if the 
   # N.B. $cpu is checked by the funtion checkFlags.
   read -e -p "What do you want to do? " choice
   case $choice in
      1|2|3|4|5|6|7|8|9|10|11|12)
         if [ -n ${timens} ]; then
            nameprod=${name1}_${timens}
            echo $nameprod 
         else
            timens="X"
            nameprod=${name1}_${timens} 
            echo $nameprod
         fi ;;
      *)
      echo " ERROR!! ERROR!! $choice not listed ERROR!! ERROR!! "
      exit
      ;;
   esac

   # Check optional dependencies [depending on the function selected]
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

   nameprod=${name1}_${timens}

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
         clusterAnalysis
         # rewrite in a "repeat" function so I can use it also with option 10
         read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
         while [ "$ramen" == "yes" ] 
         do
            repeatClusterAnalysis
         done  ;;
      9)
         pca   ;;
      10)
         clusterAnalysis && pca ;;
      11)
         sasAnalysis  ;;
      12)
         GGplot  ;;
   esac

fi
 
################ The program ends here ################
