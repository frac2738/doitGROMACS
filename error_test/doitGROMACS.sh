#!/bin/bash 

#------------------------------------------------------------------------------
#
#   File:       doitGROMACS.sh          
#   Version:    V1.6                                                    
#   Update:     08.06.15                                                  
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
#   --- doitgromacs-v.1.6
#       |--- INSTALL: useless file
#       |--- doitGROMACS_default.config: file containing "user dependable" 
#                                        variables. This file will be copied in
#                                        the working directory, where the user
#                                        is advise to check and edit it.
#       |--- Makefile: used to add the paths to the required binaries to the 
#                      configuration file.
#       |--- doitGROMACS.sh: main script
#       |--- doitGROMACS_functions  : file containing all the functions definition
#       |--- doiRGROMACS.R : R script used by the ggplot function to plot.
#       |--- doiRfunctions.R : R script containing the functions used by doiRGROMACS.R
#       |--- sas_mean.R :
#       |--- hb_mean.R :
#       |--- mdp_examples/ : directory containing examples of paramenters files
#                            used by gromacs to run simulations.
#             
#       PLANNED UPGRADES: - Comparison in R
#
#------------------------------------------------------------------------------
#
#   Usage:
#   ./doitGROMACS.sh [-u xxx] [-b xxx] [-n xxx] [-t xxx] [-k xxx] [-s xxx] 
#                    [-f xxx] [-c xxx] [-e xxx]
#   -g              - if present gromacs 5 syntax will be used
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

#---------------------------- The program begins here --------------------------

# set the directory from which the scripts are run
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
FUNCTIONS_BIN="$DIR/functions"
optionRprog="$FUNCTIONS/doitRGROMACS.R"

#today=$(date +%a-%d-%b' %T')

# source all the functions definition
source $FUNCTIONS_BIN/doitGROMACS_errorsHandling.sh
source $FUNCTIONS_BIN/doitGROMACS_flagsCheckers.sh
source $FUNCTIONS_BIN/doitGROMACS_messages.sh
source $FUNCTIONS_BIN/doitGROMACS_equilib.sh
source $FUNCTIONS_BIN/doitGROMACS_analyses.sh

while getopts "hgzb:n:t:s:f:c:e:" opt; do
 case $opt in
    h) helpMessage; exit    ;;
    z) unres=$OPTARG        ;;
    b) cpu=$OPTARG          ;;
    n) name1=$OPTARG        ;;
    t) timens=$OPTARG       ;;
    s) tpr=$OPTARG          ;;
    f) trj=$OPTARG          ;;
    c) pdb1=$OPTARG         ;;
    e) energy=$OPTARG       ;;
    g) gromacs_ver='5'      ;;
    \?) helpMessage;  exit  ;;
 esac
done
# if echo $* | grep -e "-g" -q; then

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
  # set the gromacs syntax (version 4 or 5)
  setGROMACSbinaries
  # list the options 
  doitOptions
  # check the existance of the selected option 
  read -e -p "execute option  " choice
  case $choice in
    all|emin|nvt|npt|h20|cond|rmsdfg|cluster|pca|sas|sas-sites|dssp|hb|hb-sites|ggplot|ggplot-bis|indexCreator|modvim+|meansas|meanhb)
      if [ -z ${timens} ]; then
        timens="X"
        nameprod=${name1}_${timens}
      else
        nameprod=${name1}_${timens} 
        fi ;;
    catomain)
      checkFlags_pdb; CAtoMAIN $pdb1  ;;
    split_states)
      checkFlags_pdb; split_states $pdb1  ;;
    *)
    doitOptions
    error_exit " line $LINENO, An error has occurred. Execution halted! choice '$choice' not recognised."  ;;
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
      gromCLUSTER
      # rewrite in a "repeat" function so I can use it also with option 10
      read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
      while [ "$ramen" == "yes" ] 
      do
        repeatgromCLUSTER
        i=$(($i+1))
      done  ;;
    pca)
      gromPCA   ;;
    sas)
      gromSAS  ;;
    sas-sites)
      gromSAS-sites  ;;
    dssp)
      gromDSSP ;;
    hb)
      gromHB ;;         
    hb-sites)
      gromHB-sites ;;         
    ggplot)
        GGplot ;;
    ggplot-bis) # hidden function
      sim_conditions && rmsdf && gromDSSP && GGplot ;;
    indexCreator)
        indexCreator ;;
    modvim+)
        modVim_plus  ;;
    meansas)
        mean_sas  ;;
    meanhb)
        mean_hb  ;;
  esac
  fi
#----------------------------- The program ends here ---------------------------
