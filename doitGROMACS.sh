#!/bin/bash
set -e   
#------------------------------------------------------------------------------
#
#   File:       doitGROMACS.sh          
#   Version:    V1.0.0                                                    
#   Update:     12.03.15                                                  
#
#   Copyright:  (c) Francesco Carbone, UCL, 2015
#   Author:     Francesco Carbone, UCL                                    
#   Function:   Script to execute a bunch of stuff with gromacs           
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College London
#               Gower Street
#               WC1E 6BT
#   EMail:      f.carbone.12@ucl.ac.uk
#
#------------------------------------------------------------------------------
#
#   Decription:
#   Che fa?
#
# 
#   --- doitgromacs-v.1.0
#       |--- INSTALL: useless file
#       |--- doitGROMACS_default.config: file containing "user dependable" 
#                                        variables. This file will be copied in
#                                        the working directory, where the user
#                                        is advise to check and edit.
#       |--- Makefile: used to add the paths to the required binaries to the 
#                      configuration file.
#       |--- doitGROMACS.sh: main script
#       |--- doiRGROMACS.R : R script used by the ggplot function to plot.
#       |--- mdp_examples/ : directory containing examples of paramenters files
#                            used by gromacs to run simulations.
#
#------------------------------------------------------------------------------
#
#   Usage:
#   ./doitGROMACS.sh [-u xxx] [-b xxx] [-n xxx] [-t xxx] [-k xxx] [-s xxx] 
#                    [-f xxx] [-c xxx] [-e xxx]
#   -u              - Analyse an unrest trajectory
#   -b              - Set binary location (gromacs and R)
#   -n              - Set the name
#   -t              - Set the simulation length
#   -s              - Set the tpr file
#   -f              - Set the trajectory file
#   -c              - Set the pdb file
#   -e              - Set the energy file
#
#------------------------------------------------------------------------------

#-----------
# Default error behaviour: exit with a general error message and save the output
# on doitgromacs.log. The general message may be replace while invoking 
# - error_exit "<write something>"
error_exit() {
  echo "$(date): ${1:-"
Unknown Error, execution halted "}" 2>&1 | tee -a doitgromacs.log 
  echo ""
  exit 1
}

#------------
helpMessage() {
  cat <<EOF

                        doitGROMACS.sh -  v 1.0.5  

    Copyright (c) 2013-2014, Francesco Carbone, University College London (UCL)

    This script is designed to automatise the first step of a molecular dynamics
    experiment ( solvation and equilibration) and to run some basic analyses on
    the trajectory files (.xtc), using GROMACS tools. This script allow also a
    basic analyses of UNRES (coarse grained ff) trajectories (-u).
    This script  is written with GROMACS 4.6 in mind and although it should be
    compatible with any previous versions, it is advise  to check the  validity
    of each commands before the use with an older version.
    
    At every run the scripts checks the existace of a config file located in the
    working directory (doitGROMACS.config). If a config file is not found, the 
    script uses doitGROMACS_default.config as template for a new configuration 
    file specific for the machine in which the script is run from by adding the
    correct path to both gromacs and R. If the machine is not specify with the
    "-b" flag, the config filw will contain standard paths. 
    NOTE 1: The script recognise four configurations: acrm/emerald/lappy/default.
            If gromacs or R is installed in a non standard location, the user
            have to manually edit the config file to match its system, otherwise
            doitGROMACS will not find the binaries.
     
    USAGE:  ./doitGROMACS.sh -h                     -->   HELP
            ./doitGROMACS.sh -u                     -->   UNRES analysis
            ./doitGROMACS.sh -b arg1 -n arg2 -...   -->   GROMACS analyses 


    Option   Type     Value       Description                  
    ---------------------------------------------------------------------------
    -[no]h   bool     yes         Print help info
    -u       string   txt file    Analyses of an unres trajectory
                                          
                  
                                          [ALWAYS REQUIRED in absence of -u]
    ---------------------------------------------------------------------------
    -b       string   acrm        Set the location of binaries
                                  acrm      -> Darwin building computer
                                  emerald   -> Emerald cluster
                                  lappy     -> Personal laptop
                                  default	  -> default locations
    -n       int      wt          Set the name

                                          [OPTIONAL (function dependant)]
    ---------------------------------------------------------------------------
    -t       int      200         Set the simulation length
    -s       string   .tpr        .tpr file          
    -f       string   .xtc        Trajectory file
    -c       string   .pdb        Pdb file to use to start a simulation
    -e       string   .edr        Energy file 

    NOTE 2:  In my simualtions all the output are printed in this format:
                              NAMErX_TIME
            where  NAME is the name of the mutation ( 306r ), rX is the replica
            number (r1,r2,...) and TIME is the simulation time.As a consequence
            this script takes and process output names in this form.
            
EOF
}

#------------
doitOptions() {
   cat <<EOF

                                 -----------
    ---------------------------- doitOPTIONS ----------------------------
                                 -----------

        all      - Starting from scratch ** 
        emin     - Starting from E-minimization ** 
        nvt      - Starting from NVT **
        npt      - Starting from NPT **
        h20      - Remove water from a trajectory file 
        cond     - Check the simulation conditions (U-T-P-density)
        rmsdfg   - Calculate RMSD, GYRATION RADIUS and RMSF [backbone & sidechains] 
        cluster  - Cluster analysis 
        pca      - PCA analysis 
        sas      - SAS analysis  
        dssp     - DSSP analysis
        ggplot   - Plot with ggplot (R)
        hb       - Hydrogen bonds analysis [not yet implemented]

    --------------------------------------------------------------------
                               -----------

    **  Options that require a parameter file (.mdp) that MUST be placed in the 
        same directory; the functions only accept these files:
                        
                          TEMP_min.mdp   
                          TEMP_nvt.mdp
                          TEMP_npt.mdp
                          TEMP_md.mdp

        with TEMP = temperature in Kelvin (310, 400, ...)
        The doitGROMACS tarball contains some examples of mdp file that can be 
        used after editing.

EOF
}
   
#--------------
# description : this function takes a ".stat" file from an UNRES simulation and 
#               passes it to an R script for analyses.
unresAnalyses() {
   $RSexePATH $rScripts/doitUNRES.R $unres
}

#-------
# description: Replace all the "@" with "#" using sed. Usefull to avoid conflicts
#              between R and grace.
modVim() {
  sed -i 's/@/#/g' "$1" # -i modify and save 
}

#-------
# description : It prepares all the input files required for a simulation. 
#               Using a force field and a water model, it first generates a 
#               topology file (pdb2gmx), it creates a box around the protein 
#               (editconf), solfatate (genbox) and ionised (grompp) the system.
# requirements: .pdb + .mdp 
inputs() { 
  # if -c is not set prompt and exit
  if [ -z "${pdb1}" ]; then
    helpMessage; error_exit " execution halted: a pdb file is required (-c)"
  fi 
  # create a topology file
  (echo "$optionFF"; echo "$optionWM") | $groPATH/pdb2gmx -f $pdb1             \
    -o $name1"_processed.pdb" -p topology.top -ignh -v
  # create a box
  $groPATH/editconf -f $name1"_processed.pdb" -o $name1"_inbox.pdb"           \
    -bt $optionBOX -d $optionDISTEDGE -c
  # solfatate the box
  $groPATH/genbox -cp $name1"_inbox.pdb" -cs spc216.gro -o $name1"_sol.pdb"   \
    -p topology.top
  $groPATH/grompp -f $ioniMDP -c $name1"_sol.pdb" -p topology.top             \
    -o input_ioni.tpr 
  # add ions
  read -e -p "how many ions do you need? " pioni
  $groPATH/genion -s input_ioni.tpr -p topology.top -o $name1"_ioni.pdb"      \
    -pname NA -nname CL -np $pioni
  # for negative ions use the -nn flag
}

#--------------------
# description : It minimises a structure 
# requirements: ioni.pdb + topology.top + min.mdp  
energy_minimization() {
  $groPATH/grompp -f $minMDP -c $name1"_ioni.pdb" -p topology.top              \
    -o input_min.tpr
  $groPATH/mdrun -s input_min.tpr -deffnm $name1"_min" -v
  # export the potential energy profile
  echo Potential | $groPATH/g_energy -f $name1"_min.edr" -o $name1"_potential.xvg"
  modVim $name1"_potential.xvg"
}  

#----
# description : It runs a NVT equilibration (constant volume and temperature)
# requirements: nvt.mdp + min.gro + topology.top
nvt() {
   $groPATH/grompp -f $nvtMDP -c $name1"_min.gro" -p topology.top              \
      -o input_nvt.tpr
   $groPATH/mdrun -s input_nvt.tpr -deffnm $name1"_nvt" -v
   # export the temperature profile 
   echo Temperature | $groPATH/g_energy -f $name1"_nvt.edr" -o $name1"_temperature.xvg"
   modVim $name1"_temperature.xvg"
} 

#----
# description : It runs a NPT equilibration (constant pressure)
# requirements: npt.mdp + nvt.gro + topology.top
npt() {
   $groPATH/grompp -f $nptMDP -c $name1"_nvt.gro" -p topology.top              \
      -o input_npt.tpr
   $groPATH/mdrun -s input_npt.tpr -deffnm $name1"_npt" -v  
   # export the pressure profile 
   echo Pressure | $groPATH/g_energy -f $name1"_npt.edr" -o $name1"_pressure.xvg"
   modVim $name1"_pressure.xvg"
   # export the density profile
   echo Density | $groPATH/g_energy -f $name1"_npt.edr" -o $name1"_density.xvg"
   modVim $name1"_density.xvg"
} 

#---------------
# description : Given an energy file, it extracts energy, temerature, pressure
#               and density profiles.
# requirements: energy.edr
sim_conditions() {
   # if -e is not set prompt and exit
   if [ -z "${energy}" ]; then
      helpMessage; error_exit " execution halted: an energy file is required (-e)"
   fi

   # export potential energy, temperature, pressure and density profiles 
   echo Potential | $groPATH/g_energy -f $energy -o $nameprod"_potential.xvg"
   modVim $nameprod"_potential.xvg"
   echo Temperature | $groPATH/g_energy -f $energy -o $nameprod"_temperature.xvg"
   modVim $nameprod"_temperature.xvg"
   echo Pressure | $groPATH/g_energy -f $energy -o $nameprod"_pressure.xvg"
   modVim $nameprod"_pressure.xvg"
   echo Density | $groPATH/g_energy -f $energy -o $nameprod"_density.xvg"
   modVim $nameprod"_density.xvg"
   # to plot using ggplot check option the function "GGplot"
}

checkTflags() {
  if [[ -z "${tpr+x}" || -z "${trj+x}" ]]; then
    helpMessage
    error_exit " execution halted: tpr and trj files are required (-s -f)"
  fi
}

#----------
# description : It removes the water molecules from a trajectory (.xtc) and 
#               removes the pbc effects (pbc = periodic buonday conditions).
# requirements: .tpr + .trj 
clean_trj() {
   # if -s and -f are not set prompt and exit
   checkTflags
   # removing water molecules from the trajectory file
   echo Protein | $groPATH/trjconv -s $tpr -f $trj                             \
      -o $nameprod"_only.xtc"
   # removing water molecules from the tpr file 
   echo Protein | $groPATH/tpbconv -s $tpr -o $nameprod"_only.tpr"
   # creating a protein only .gro file 
   echo Protein | $groPATH/trjconv -s $nameprod"_only.tpr"                     \
      -f $nameprod"_only.xtc" -o $nameprod"_only.gro" -dump 0
   # account for the periodicity (nojump)
   echo Protein | $groPATH/trjconv -s $nameprod"_only.tpr"                     \
      -f $nameprod"_only.xtc" -o $nameprod"_nojump.xtc" -pbc nojump
   # account for the periodicity (fitting) 
   (echo "Backbone"; echo "Protein") | $groPATH/trjconv -s $nameprod"_only.tpr"\
      -f $nameprod"_nojump.xtc" -o $nameprod"_fit.xtc" -fit rot+trans
   # remove intermediate files and rename them in a less complicated way
   rm $nameprod"_nojump.xtc" $nameprod".xtc" $nameprod"_only.xtc"
   mv $nameprod"_fit.xtc" $nameprod".xtc"
   mv $nameprod"_only.tpr" $nameprod".tpr"
   mv $nameprod"_only.gro" $nameprod".gro"
}

#------
# description : It calculates rmsd (backbone), gyration radius and rmsf (both 
#               for the backbone and the sidechains) 
# requirements: .tpr + .xtc 
rmsdf() {
   checkTflags
   # calculating the RMSD 
   (echo "$optionRMSD"; echo "$optionRMSD") | $groPATH/g_rms -s $tpr -f $trj   \
      -o $nameprod"_rmsd.xvg" -tu ns
   modVim $nameprod"_rmsd.xvg"
   # calculating the radius of Gyration 
   echo $optionGYRATION | $groPATH/g_gyrate -s $tpr -f $trj                    \
      -o $nameprod"_rgyration.xvg"
   modVim $nameprod"_rgyration.xvg"
   echo $optionRMSFb | $groPATH/g_rmsf -s $tpr -f $trj                         \
      -o $nameprod"_rmsf_bb.xvg" -oq $nameprod"_rmsf_bb.pdb" -res
   modVim $nameprod"_rmsf_bb.xvg" 
      # res=averages for each residues
   echo $optionRMSFsc | $groPATH/g_rmsf -s $tpr -f $trj                        \
      -o $nameprod"_rmsf_sc.xvg" -oq $nameprod"_rmsf_sc.pdb" -res 
   modVim $nameprod"_rmsf_sc.xvg" 
}

#----------------
# description : read the name of the function
# requirements: .tpr + .xtc
clusterAnalysis() {
   # if -s and -f are not set prompt and exit
   checkTflags
   if [ ! -d ./clusters_$name1 ] ; then
      mkdir clusters_$name1
   fi
   cd clusters_$name1
   # check if the rmsd matrix exists
   while [ ! -f ./rmsd-matrix.xpm ]
   do
      # create the rmsd-matrix (using the backbone for the calculation)
      (echo "$optionCLUSTER"; echo "$optionCLUSTER") | $groPATH/g_rms -s ../$tpr -f ../$trj   \
         -f2 ../$trj -m rmsd-matrix.xpm -dist rmsd-distribution.xvg -tu ns     \
         -skip $optionSKIPcluster -b $optionSTARTimeNS
      # improve rmsd matrix plot
      $groPATH/xpm2ps -f rmsd-matrix.xpm -o rmsd-matrix.eps
   done
   #check the distribution file and decide the cutoff
   xmgrace rmsd-distribution.xvg
   read -e -p "Which cutoff do you want to use? " cutoff
   # cluster analysis on the rmsd matrix (Using the backbone for the calculation)
   (echo "$optionCLUSTER"; echo "$optionCLUSTER") | $groPATH/g_cluster -s ../$tpr -f ../$trj  \
      -dm rmsd-matrix.xpm -o clusters.xpm -sz clusters-size.xvg                \
      -clid clusters-ovt.xvg -cl clusters.pdb -cutoff $cutoff                  \
      -method $optionCLUSTERMETHOD -tu ns -skip $optionSKIPcluster -b $optionSTARTimeNS
   # to visualize in pymol use
   # split_states clusters
   # delete clusters
   # dss
   # show cartoon
   cd ..
}

#----------------------
# description : The name says everything
# requirements: see "clusterAnalysis"
repeatClusterAnalysis() {
   mv clusters_$name1 clusters"$i"_$name1
   clusterAnalysis
   i++
   read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
}

#----
# description : The name says everything
# requirements: .tpr + .xtc
pca() {
  # if -s and -f are not set prompt and exit
  checkTflags
  if [ ! -d ./PCA_$name1 ] ; then
     mkdir PCA_$name1
  fi
  cd PCA_$name1
  # check if the covariance matrix exists
  while [ ! -f ./covariance.xpm ]
  do
     # calculating the covariance matrix (C-alpha), eigenvalues and eigenvectors
     (echo "$optionPCA"; echo "$optionPCA") | $groPATH/g_covar -s ../$tpr -f ../$trj   \
        -o eigenvalues.xvg -v eigenvectors.trr -xpma covariance.xpm -tu ns    \
        -dt $optionDTpca -b $optionSTARTimeNS
     # improve covariance matrix ploting
     $groPATH/xpm2ps -f covariance.xpm -o covariance.eps
  done
  xmgrace eigenvalues.xvg
  # determinare su quali autovalori lavorare
  read -e -p "how many eigenvalues do you want to use for the analysis? " range
  (echo "$optionPCA"; echo "$optionPCA") | $groPATH/g_anaeig -v eigenvectors.trr       \
     -s ../$tpr -f ../$trj -first 1 -last $range -proj "projection-1"$range".xvg" \
     -tu ns -b $optionSTARTimeNS
  for ((i = 1; i <= range; i=i+1))
  do
     (echo "$optionPCA"; echo "$optionPCA") | $groPATH/g_anaeig -v eigenvectors.trr \
        -s ../$tpr -f ../$trj -first $i -last $i -nframes 100                 \
        -extr "ev"$i".pdb" -filt "ev"$i".xtc" -proj "ev"$i".xvg" -tu ns       \
        -b $optionSTARTimeNS
     $groPATH/g_analyze -f "ev"$i".xvg" -cc "ev"$i"-cc.xvg" -ac "ev"$i"-ac.xvg"
  done
  # calculate 2d projections and FES
  (echo "$optionPCA"; echo "$optionPCA") | $groPATH/g_anaeig -v eigenvectors.trr   \
     -s ../$tpr -f ../$trj -first 1 -last 2 -2d 2d-12.xvg
  $groPATH/g_sham -f 2d-12.xvg -ls gibbs-12.xpm -notime
  $groPATH/xpm2ps -f gibbs-12.xpm -o gibbs-12.eps -rainbow red
  perl /export/francesco/Dropbox/scripts/Rscripts/Rplots/heatmap.pl gibbs-12.xpm
  $RscriptEXE $optionRprog -hm=pes_profile.txt  
  mv pes.png ../$name1"_pes.png"
  cd ..
}

#------------
# description : The name says everything
# requirements: .tpr + .xtc
sasAnalysis() {
   # if -s and -f are not set prompt and exit
   checkTflags
   # create the directory
   if [ ! -d ./sas_$name1 ] ; then
      mkdir sas_$name1
   fi
   cd sas_$name1
   (echo "$optionSAS"; echo "$optionSAS") | $groPATH/g_sas -s ../$tpr -f ../$trj \
      -o $name1"_area.xvg" -tv $name1"_volume.xvg" -or $name1"_resarea.xvg"      \
      -dt $optionDTsas -b $optionSTARTimePS -probe $optionPROBE
   # probe 0.7 nm perchè 16A -> 1.6 nm -> raggio 0.7 nm
   cd ..
}

#-------
# description :
# requirements: .tpr + .xtc
gromHB() {
  checkTflags
  (echo "optionHB; echo "optionHB"")|$groPATH/g_hbond -s $trp -f $trj         \
    -num $nameprod"_hb_count.xvg" -dist $nameprod"_hb_dist.xvg" -b $optionSTARTimeNS
    -hbm $nameprod"_hb_matrix" -tu ns -dt $optionDT -b $optionSTARTimeNS 
}

#---------
# description :
# requirements: .tpr + .xtc
gromDSSP() {
  # if -s and -f are not set prompt and exit
  checkTflags
  echo $optionDSSP | $groPATH/do_dssp -ver 1 -f $trj -s $tpr -o $nameprod"_ss.xpm" \
    -sc $nameprod"_ss_count.xvg" -tu ns -dt $optionDTdssp -b $optionSTARTimeNS  
  modVim $nameprod"_ss_count.xvg"
}

#-------
GGplot() {
  # function that calls the R script doitRGROMACS.R to plot in ggplot: 
  # rmsd - gyration radius - rmsf - simulation conditions - structure analysis
  if [[ -f $optionRprog ]]; then
    $RscriptEXE $optionRprog -d=$nameprod"_rmsd.xvg"                           \
      -g=$nameprod"_rgyration.xvg" -ss=$nameprod"_ss_count.xvg"                \
      -x=$nameprod"_density.xvg" -t=$nameprod"_temperature.xvg"                \
      -u=$nameprod"_potential.xvg" -p=$nameprod"_pressure.xvg"                 \
      -fb=$nameprod"_rmsf_bb.xvg" -fsc=$nameprod"_rmsf_sc.xvg"
    # rename the outputs
  else 
    error_exit " the function "GGplot" requires a R script located in $DIR."
  fi 
}

checkFlags() {
  # check the definition of -b and -n 
  if [[ -z "${cpu+x}" || -z "${name1+x}" ]]; then
    helpMessage
    error_exit "Execution halted! the script was called without any flags." 
  fi 
}

#-------------------------------------------------------------------------------
#---------------------------- The program begins here --------------------------

# set the directory from which the script is run & the R script
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
optionRprog="$DIR/doitRGROMACS.R"

#today=$(date +%a-%d-%b' %T')

while getopts "hu:b:n:t:s:f:c:e:" opt; do
 case $opt in
    h) helpMessage; exit    ;;
    u) unres=$OPTARG        ;;
    b) cpu=$OPTARG          ;;
    n) name1=$OPTARG        ;;
    t) timens=$OPTARG       ;;
    s) tpr=$OPTARG          ;;
    f) trj=$OPTARG          ;;
    c) pdb1=$OPTARG         ;;
    e) energy=$OPTARG       ;;
    \?) helpMessage;  exit  ;;
 esac
done

# if -u exists call unres(), otherwise execute gromacs
if [ -n "${unres}" ]; then
  unresAnalyses; exit;   
else
  checkFlags

  # check existance of CONFIG_FILE and source it or create a new one
  export CONFIG_FILE="$(find . -maxdepth 1 -name doitGROMACS.config)"
  # Source configuration file
  if [[ -f $CONFIG_FILE ]]; then
    . $CONFIG_FILE ;
 else
    echo "Configuration file not found, a new file will be created";
    cp $DIR/doitGROMACS_default.config ./doitGROMACS.config
    case $cpu in
      acrm | emerald | lappy) 
      make -f $DIR/Makefile $cpu
      . doitGROMACS.config
      echo "
------------------------------ executables found ----------------------------
executables locations specific for $cpu have been written on doitGROMACS.conf
-----------------------------------------------------------------------------
      "	;;
      *)	
      make -f $DIR/Makefile standard
      . doitGROMACS.config
      echo "
---------------------- no executables found ------------------------
standard executables locations have been written on doitGROMACS.conf
--------------------------------------------------------------------
      "	;;
    esac 
  fi
  # list the options 
  doitOptions

  # check the existance of the selected option 
  read -e -p "execute option  " choice
  case $choice in
    all|emin|nvt|npt|h20|cond|rmsdfg|cluster|pca|sas|dssp|hb|ggplot|ggplot-bis)
      if [ -z ${timens} ]; then
        timens="X"
        nameprod=${name1}_${timens}
      else
        nameprod=${name1}_${timens} 
        fi ;;
    *)
    #echo " ERROR!! ERROR!! $choice not listed ERROR!! ERROR!! " >> doitGROMACS.log 2>&1
    doitOptions
    error_exit " Execution halted! choice '$choice' not recognised."  ;;
  esac

  case $choice in
    all )
      inputs && energy_minimization && nvt && npt  ;;
    emin)
      energy_minimization && nvt && npt   ;;
    nvt)
      nvt && npt  ;;
    npt)
      npt   ;;
    h20)
      clean_trj   ;;
    cond) 
      sim_conditions ;;
    rmsdfg)
      rmsdf ;;
    cluster)
      i=1
      clusterAnalysis
      # rewrite in a "repeat" function so I can use it also with option 10
      read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
      while [ "$ramen" == "yes" ] 
      do
        repeatClusterAnalysis
      done  ;;
    pca)
      pca   ;;
    sas)
      sasAnalysis  ;;
    dssp)
      gromDSSP ;;
    hb)
      gromHB ;;         
    ggplot)
        GGplot ;;
    ggplot-bis) # hidden function
      sim_conditions && rmsdf && gromDSSP && GGplot
      mv rmsd.png $name1"_rmsd.png"
      mv rgyr.png $name1"_rgyr.png" 
      mv rmsf.png $name1"_rmsf.png"
      mv simcond.png $name1"_simcond.png"
      mv ss.png $name1"_ss.png" ;;
  esac
  fi
#----------------------------- The program ends here ---------------------------
