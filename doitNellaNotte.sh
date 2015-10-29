#!/bin/bash
set -e

# 450
cd /media/g6pd8/simulations_data/450K/trj_450

echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n wt -t 500 -s wt_500.tpr -f wt_500.xtc -e wt_500.edr

echo "sas" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n Ameno -t 500 -s Ameno_500.tpr -f Ameno_500.xtc -e Ameno_500.edr
echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n Ameno -t 500 -s Ameno_500.tpr -f Ameno_500.xtc -e Ameno_500.edr

echo "sas" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n 306r -t 500 -s 306r_500.tpr -f 306r_500.xtc -e 306r_500.edr
echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n 306r -t 500 -s 306r_500.tpr -f 306r_500.xtc -e 306r_500.edr

echo "sas" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n 338e -t 500 -s 338e_500.tpr -f 338e_500.xtc -e 338e_500.edr
echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n 338e -t 500 -s 338e_500.tpr -f 338e_500.xtc -e 338e_500.edr

# 470
cd /media/g6pd8/simulations_data/470K/trj_470

echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n wt -t 500 -s wt_500.tpr -f wt_500.xtc -e wt_500.edr
echo "sas" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n wt -t 500 -s wt_500.tpr -f wt_500.xtc -e wt_500.edr

echo "sas" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n Ameno -t 500 -s Ameno_500.tpr -f Ameno_500.xtc -e Ameno_500.edr
echo "sas-sites" | /export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/doitGROMACS.sh -g -b acrm -n Ameno -t 500 -s Ameno_500.tpr -f Ameno_500.xtc -e Ameno_500.edr

