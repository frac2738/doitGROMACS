all: no_option
	@echo "Specify machine [acrm-emerald-lappy-default]"

no_option:
   
acrm: 
	#cp doitGROMACS_default.config doitGROMACS.config
	echo "# Gromacs"	>> doitGROMACS.config
	echo "groPATH='/acrm/usr/local/apps/gromacs/bin'"	>> doitGROMACS.config
	echo "# R"														>> doitGROMACS.config
	echo "REXE='/export/francesco/R-3.1.0/bin/R'"     	>> doitGROMACS.config	
	echo "RscriptEXE='/export/francesco/R-3.1.0/bin/Rscript'" >> doitGROMACS.config
	echo "# R scripts directories"							>> doitGROMACS.config
	echo "rScriptsDIR='/export/francesco/Dropbox/scripts/Rscripts/Rplots'" >> doitGROMACS.config

emerald: 
	#cp doitGROMACS_default.config doitGROMACS.config
	echo "# Gromacs"										>> doitGROMACS.config
	echo "groPATH='/apps/gromacs/4.6.3/bin'"		>> doitGROMACS.config

lappy: 
	#cp doitGROMACS_default.config doitGROMACS.config
	echo "# Gromacs"										>> doitGROMACS.config
	echo "groPATH='/usr/local/gromacs/bin'"  		>> doitGROMACS.config
	echo "# R"												>> doitGROMACS.config
	echo "REXE='/usr/bin/R'"     						>> doitGROMACS.config                            
	echo "RscriptEXE='/usr/bin/Rscript'"			>> doitGROMACS.config
	echo "# R scripts directories"					>> doitGROMACS.config
	echo "rScriptsDIR='~/Dropbox/scripts/Rscripts/Rplots'"	>> doitGROMACS.config

standard:
	#cp doitGROMACS_default.config doitGROMACS.config
	echo "# Gromacs"										>> doitGROMACS.config
	echo "groPATH='/usr/local/gromacs/bin'"		>> doitGROMACS.config
	echo "# R"												>> doitGROMACS.config	
	echo "REXE='/usr/bin/R' "							>> doitGROMACS.config                               
	echo "RscriptEXE='/usr/bin/Rscript'" 			>> doitGROMACS.config
	echo "# R scripts directories"					>> doitGROMACS.config
	echo "rScriptsDIR='./'"								>> doitGROMACS.config

clean:
	echo "Copyright (c) 2013-2014, Francesco Carbone, University College London (UCL)"
