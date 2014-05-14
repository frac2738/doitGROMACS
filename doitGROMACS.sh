#!/bin/bash
set -e   
###############################################################################
#     Version:      V1.0.3                                                    #
#     Last update:  19.02.14                                                  #
#     Author:       Francesco Carbone                                         #
#     Description:  Script to execute a bunch of stuff with gromacs           #
#     Updates :     - replaced "if" with "case" statement                     #
#                   - added the rmsdf function                                #
#                   - added the "sim_conditions" function                     #
#                   - added the "patches" function [testing]                  #
###############################################################################

# functions definition

function inputs {
   ls
   read -e -p "which pdb should be used as input? " pdb1

   echo "-------- TOPOLOGY CREATION --------"
   $path/pdb2gmx -f $pdb1 -o $name1"_processed.pdb" -p topology.top -ignh -v

   echo "-------- BOX DEFINITION --------"
   read -e -p "select the type of box? " boxtype
   read -e -p "select the protein-box distance? (nm) " distedge
   $path/editconf -f $name1"_processed.pdb" -o $name1"_inbox.pdb" -bt $boxtype \
      -d $distedge -c

   echo "-------- SOLVATION --------"
   $path/genbox -cp $name1"_inbox.pdb" -cs spc216.gro -o $name1"_sol.pdb"      \
      -p topology.top
   $path/grompp -f G6PD_ionization.mdp -c $name1"_sol.pdb" -p topology.top     \
      -o input_ioni.tpr
   read -e -p "how many ions do you need? " pioni
   $path/genion -s input_ioni.tpr -p topology.top -o $name1"_ioni.pdb"         \
      -pname NA -nname CL -np $pioni
   # for negative ions use the -np flag
}

function energy_minimization {
   echo "-------- ENERGY MINIMIZATION --------"
   $path/grompp -f $temp1"_min.mdp" -c $name1"_ioni.pdb" -p topology.top       \
      -o input_min.tpr
   $path/mdrun -s input_min.tpr -deffnm $name1"_min" -v  
   # for X processor use:
   # $path/mpirun -np X mdrun_mpi -s input_min.tpr -deffnm $name1"_min" -v
   # plot the potential energy profile
   echo Potential | $path/g_energy -f $name1"_min.edr" -o $name1"_potential.xvg"
}  

function nvt {
   echo "-------- NVT --------"
   $path/grompp -f $temp1"_nvt.mdp" -c $name1"_min.gro" -p topology.top        \
      -o input_nvt.tpr
   $path/mdrun -s input_nvt.tpr -deffnm $name1"_nvt" -v
   # plot the temperature profile 
   echo Temperature | $path/g_energy -f $name1"_nvt.edr" -o $name1"_temperature.xvg"
} 

function npt {
   echo "-------- NPT --------"
   $path/grompp -f $temp1"_npt.mdp" -c $name1"_nvt.gro" -p topology.top        \
      -o input_npt.tpr
   $path/mdrun -s input_npt.tpr -deffnm $name1"_npt" -v  
   # plot the pressure profile 
   echo Pressure | $path/g_energy -f $name1"_npt.edr" -o $name1"_pressure.xvg"
   # plot the density profile
   echo Density | $path/g_energy -f $name1"_npt.edr" -o $name1"_density.xvg"
} 

function sim_conditions {
   # calculate system potential energy 
   echo Potential | $path/g_energy -f $nameprod".edr" -o $nameprod"_potential.xvg"
   # calculate system temperature 
   echo Temperature | $path/g_energy -f $nameprod".edr" -o $nameprod"_temperature.xvg"
   # calculate system pressure 
   echo Pressure | $path/g_energy -f $nameprod".edr" -o $nameprod"_pressure.xvg"
   # calculate system density
   echo Density | $path/g_energy -f $nameprod".edr" -o $nameprod"_density.xvg"
}

function clean_trj {
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
   # remove intermediate files
   rm $nameprod"_nojump.xtc" $nameprod".xtc" $nameprod"_only.xtc"
   mv $nameprod"_fit.xtc" $nameprod".xtc"
   mv $nameprod"_only.tpr" $nameprod".tpr"
   mv $nameprod"_only.gro" $nameprod".gro"
}

function rmsdf {
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

function cluster_analysis {
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
   read -e -p "Which method do you want to use? [gromos - linkage] " method
   # cluster analysis on the rmsd matrix (Using the backbone for the calculation)
   (echo "Backbone"; echo "Backbone") | $path/g_cluster -s ../$tpr -f ../$trj  \
      -dm rmsd-matrix.xpm -o clusters.xpm -sz clusters-size.xvg                \
      -clid clusters-ovt.xvg -cl clusters.pdb -cutoff $cutoff -method $method  \
      -tu ns -skip $skip
   # repeat the analysis using a different -method to improve results
   #
   # to visualize in pymol use
   # split_states clusters
   # delete clusters
   # dss
   # show cartoon
   cd ..
}

function pca {
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

function patches {
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
   /home/bsm3/zcbtfo4/public/automatic_patches -r 8 -t normal -l $nameprod"_patchlog" -o $nameprod"_patches" -
   cd ..
}

function plot {
   write something
}

################ THE PROGRAM BEGINS HERE ################

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
echo " 10 - Do 8 - 9 "
echo " 11 - Patch analysis [soon] "
echo " 12 - Plot [soon] "
echo "-----------------------------------------------------------"

read -e -p "What do you want to do? " choice
case $choice in
   1|2|3|4|5|6|7|8|9|10|11|12)
   ;; 
   *)
   echo "ERROR!! ERROR!! ERROR!! $choice not listed" ERROR!! ERROR!! ERROR!!
   exit
   ;;
esac

# define the path to gromacs
echo " acrm17  => /acrm/usr/local/apps/gromacs/bin"
echo " Emerald => /apps/gromacs/4.6.3/bin"
echo " bear    => /usr/local/gromacs/bin"
read -e -p "Enter the path to gromacs: " path   # -i "/acrm/usr/local/apps/gromacs/bin" path

echo "All the output names will be written following this rules: NAME_rX_TIME"
read -e -p "Insert NAME " name1
read -e -p "Insert TIME (ns) " timens
read -e -p "Insert REPLICA NUMBER " repli1
nameprod="${name1}_${repli1}_${timens}"
nameprod2="${name1}_${repli1}"

if [ $choice = 1 -o $choice = 2 -o $choice = 3 -o $choice = 4 ]; then 
   read -e -p "What temperature are you using? " temp1
fi

if [ $choice = 8 -o $choice = 9 -o $choice = 10 -o $choice = 11 ]; then
   ls
   read -e -p "Select the .tpr file : " tpr
   read -e -p "Select the trajectory file: " trj
fi

case $choice in
   1 )
      inputs && energy_minimization && nvt && npt # ;; x separare comandi/ && se il 2° MUST wait prima di essere eseguito
      ;;             # with ";&" the next block is executer without testing
                     # with ";;&" il programma passa al blocco successivo ma lo testa anche se c'è già stato il match 
   2)
      energy_minimization && nvt && npt
      ;;
   3)
      nvt && npt
      ;;
   4)
      npt
      ;;
   5)
      clean_trj
      ;;
   6) 
      sim_conditions
      ;;
   7)
      rmsdf
      ;;
   8)
      i=1
      cluster_analysis
      read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
      while [ "$ramen" == "yes" ] 
      do
         mv clusters_$nameprod2 clusters"$i"_$nameprod2
         cluster_analysis
         read -e -p "Do you want to rerun the analysis with a different method? [yes/no] " ramen
         i++
      done
      ;;
   9)
      pca
      ;;
   10)
      cluster_analysis && pca
      ;;
   11)
      patches
      ;;
   12)
      plot
      ;;
esac
 
################ THE PROGRAM ENDS HERE ################
