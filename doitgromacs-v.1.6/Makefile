all: no_option

no_option:
	@echo "Specify machine [acrm-emerald-lappy-standard]"
   
acrm: 
	echo "# Gromacs"		>> doitGROMACS.config
	echo "groPATH4='/acrm/usr/local/apps/gromacs/bin'"	>> doitGROMACS.config	
	echo "groPATH5='/acrm/usr/local/apps/gromacs-5.0.4/bin'"	>> doitGROMACS.config	
	echo "# R for acrm17"	>> doitGROMACS.config
	echo "REXE='/export/francesco/R-3.1.0/bin/R'"	>> doitGROMACS.config	
	echo "RscriptEXE='/export/francesco/R-3.1.0/bin/Rscript'"	>> doitGROMACS.config
	echo ""	>> doitGROMACS.config
	echo "# R for any other machines"	>> doitGROMACS.config
	echo "#REXE='/usr/bin/R'"	>> doitGROMACS.config	
	echo "#RscriptEXE='/usr/Rscript'"	>> doitGROMACS.config
  
emerald: 
	echo "# Gromacs"	>> doitGROMACS.config
	echo "#groPATH4='/apps/gromacs/4.6.3/bin'"	>> doitGROMACS.config
	echo "groPATH4='/apps/gromacs/4.6.7/bin'"	>> doitGROMACS.config
	echo "groPATH5='/apps/gromacs/5.0.4/bin'"	>> doitGROMACS.config

lappy: 
	echo "# Gromacs"	>> doitGROMACS.config
	echo "groPATH4='/usr/local/gromacs/bin'"	>> doitGROMACS.config
	echo "# R"	>> doitGROMACS.config	
	echo "REXE='/usr/bin/R'"	>> doitGROMACS.config          
	echo "RscriptEXE='/usr/bin/Rscript'"	>> doitGROMACS.config

standard:
	echo "# Gromacs"	>> doitGROMACS.config
	echo "groPATH4='/usr/local/gromacs/bin'"	>> doitGROMACS.config
	echo "# R"	>> doitGROMACS.config	
	echo "REXE='/usr/bin/R' "	>> doitGROMACS.config         
	echo "RscriptEXE='/usr/bin/Rscript'"	>> doitGROMACS.config

clean:
	@echo ""
	@echo "	Copyright (c) 2013-2015, Francesco Carbone, University College London (UCL)"
	@echo ""
