#!/bin/bash
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

modVim_plus() {
  read -e -p "che file vuoi? " file1 
  modVim $file1
}

inputs() { 
  # if -c is not set prompt and exit
  checkFlags_pdb 
  # create a topology file
  (echo "$optionFF"; echo "$optionWM") | $groPATH/$pdb2gmx -f $pdb1            \
    -o $name1"_processed.pdb" -p topology.top -ignh -v || checkExitCode 
  # create a box
  $groPATH/$editconf -f $name1"_processed.pdb" -o $name1"_inbox.pdb"           \
    -bt $optionBOX -d $optionDISTEDGE -c || checkExitCode 
  # solfatate the box
  $groPATH/$genbox -cp $name1"_inbox.pdb" -cs spc216.gro -o $name1"_sol.pdb"   \
    -p topology.top || checkExitCode
  $groPATH/$grompp -f $ioniMDP -c $name1"_sol.pdb" -p topology.top             \
    -o input_ioni.tpr || checkExitCode 
  read -e -p "Specify the number of ions to be add added to the system in the form of [+/- n°]
(e.g. + 12 or - 23) " number_ioni
  split_ioni=( $number_ioni )
  if [ ${split_ioni[0]} == "+" ]; then
    optionIONSpos=${split_ioni[1]}
    optionIONSneg='0'
  else 
    optionIONSneg=${split_ioni[1]}
    optionIONSpos='0'
  fi 
  # add ions
  $groPATH/$genion -s input_ioni.tpr -p topology.top -o $name1"_ioni.pdb"      \
    -pname NA  -nname CL -np $optionIONSpos -nn $optionIONSneg ; checkExitCode
  # [-np] for positive and [-nn] for negative 
} > >(tee doitgromacs_inputs.log) 2> >(tee doitgromacs_inputs.err >&2)

energy_minimization() {
  $groPATH/$grompp -f $minMDP -c $name1"_ioni.pdb" -p topology.top              \
    -o input_min.tpr
  $groPATH/mdrun -s input_min.tpr -deffnm $name1"_min" -v
  # export the potential energy profile
  echo Potential | $groPATH/$g_energy -f $name1"_min.edr" -o $name1"_potential.xvg"
  modVim $name1"_potential.xvg"
} &> >(tee doitgromacs_emin.log) >&2

nvt() {
   $groPATH/$grompp -f $nvtMDP -c $name1"_min.gro" -p topology.top              \
      -o input_nvt.tpr
   $groPATH/mdrun -s input_nvt.tpr -deffnm $name1"_nvt" -v
   # export the temperature profile 
   echo Temperature | $groPATH/$g_energy -f $name1"_nvt.edr" -o $name1"_temperature.xvg"
   modVim $name1"_temperature.xvg"
} &> >(tee doitgromacs_nvt.log) >&2

npt() {
   $groPATH/$grompp -f $nptMDP -c $name1"_nvt.gro" -p topology.top              \
      -o input_npt.tpr -maxwarn 1
   $groPATH/mdrun -s input_npt.tpr -deffnm $name1"_npt" -v  
   # export the pressure profile 
   echo Pressure | $groPATH/$g_energy -f $name1"_npt.edr" -o $name1"_pressure.xvg"
   modVim $name1"_pressure.xvg"
   # export the density profile
   echo Density | $groPATH/$g_energy -f $name1"_npt.edr" -o $name1"_density.xvg"
   modVim $name1"_density.xvg"
} &> >(tee doitgromacs_npt.log) >&2

