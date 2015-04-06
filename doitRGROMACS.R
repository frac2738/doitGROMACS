#! /export/francesco/R-3.1.0/bin/Rscript
#------------------------------------------------------------------------------
#
# File:       Rplot_gromacs.R          
# Version:    V1.0.0                                                    
# Update:     26.02.2015                                                  
#
# Copyright:  (c) Francesco Carbone, UCL, 2015
# Author:     Francesco Carbone, UCL                                    
# Function:   Script to plot xvg files           
# Address:    Institute of Structural and Molecular Biology
#             Division of Biosciences
#             University College London
#             Gower Street
#             WC1E 6BT
# EMail:      f.carbone.12@ucl.ac.uk
#
#------------------------------------------------------------------------------
#
# Decription:
# This script pplots in ggplot
#
# NOTE 1: This script is design to work with xvg derivated from a MD experiment
#         on G6PD enzyme, which means that:
#         - the system has two chains
#         - the system has 958 residues
#         - all the xvg files have been generated using GROMACS in grace format
#
#         If a system differ for all, of some of the conditions listed above,
#         wait...I'll release a better version at some point in the future.
#
# Known bugs: a) ggplot fonts don't render correctly in machines other than acrm
#
#------------------------------------------------------------------------------
#
# Usage:
# Rscript doitRGROMACS.R [-x=xx] [-t=xx] [-u=xx] [-p=xx] [-d=xx] [-g=xx]
#                  [-ss=xx] [-fb=xx] [-fsc=xx] [-help]
#
# -help   -  ask for help
# -h      -  read hight
# -w      -  read width
# -dpi    -  read dpi
# -x      -  read density file (.xvg)
# -t      -  read temperature file (.xvg)
# -u      -  read potential energy file (.xvg)
# -p      -  read pressure file (.xvg)
# -d      -  read rmsd file (.xvg)
# -g      -  read gyration radius file (.xvg)
# -ss     -  read secondary structure file (.xvg)
# -fb     -  read rmsf file [backbone] (.xvg)
# -fsc    -  read rmsf file [side chains] (.xvg)
# -hm     -  read pes projection (heat map)
#------------------------------------------------------------------------------

#---------------- Arguments parsing ----------------
# Collect the arguments
args <- commandArgs(TRUE)

# call the help when no arguments are passed
if(length(args)<1) {
   args <- c("-help")
}
print(args)    # check the args vector

## Help section
if ("-help" %in% args) {
  cat("
    Usage:   Rscript doitRGROMACS.R [-x=xx] [-t=xx] [-u=xx] [-p=xx] [-d=xx]
                  [-g=xx] [-ss=xx] [-fb=xx] [-fsc=xx] [-help]

      -help   - ask for help
      -o      - desired output
      -u      - read potential energy file (.xvg)
      -t      - read temperature file (.xvg)
      -p      - read pressure file (.xvg)
      -x      - read density file (.xvg)
      -d      - read rmsd file (.xvg)
      -g      - read gyration radius file (.xvg)
      -ss     - read secondary structure file (.xvg)
      -fb     - read rmsf file [backbone] (.xvg)
      -fsc    - read rmsf file [side chains] (.xvg)
      -hm     - read pes projection (heat map)
      -h      - set image height in
      -w      - set image width in
      -dpi    - set image dpi

      NOTE1 :  flags x,t,u,p are read together, the absence of one (or more) of
               them will cause the function to exit without plotting anything. 
      
")
   q(save="no")
}

## Parse the arguments in the form of -arg=value
parseArgs <- function(x) strsplit (sub("^-","",x), "=")   # split flags from values
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
#argsL <- as.list(as.character(argsDF$V2))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- (argsDF$V1)

## default behaviours
if(is.null(argsL$x)) print("density file not supplied") else file.density <- argsDF$V2[which(argsDF$V1 == "x")] 
if(is.null(argsL$t)) print("temperature file not supplied") else file.temperature <- argsDF$V2[which(argsDF$V1 == "t")] 
if(is.null(argsL$u)) print("energy file not supplied") else file.potential <- argsDF$V2[which(argsDF$V1 == "u")] 
if(is.null(argsL$p)) print("pressure file not supplied") else file.pressure <- argsDF$V2[which(argsDF$V1 == "p")] 
if(is.null(argsL$d)) print("rmsd file not supplied") else file.rmsd <- argsDF$V2[which(argsDF$V1 == "d")]
if(is.null(argsL$g)) print("rgyration file not supplied") else file.rgyr <- argsDF$V2[which(argsDF$V1 == "g")]
if(is.null(argsL$ss)) print("secondary structure count file not supplied") else file.ss <- argsDF$V2[which(argsDF$V1 == "ss")]
if(is.null(argsL$fb)) print("rmsf bb file not supplied") else file.rmsf.bb <- argsDF$V2[which(argsDF$V1 == "fb")]
if(is.null(argsL$fsc)) print("rmsf sc file not supplied") else file.rmsf.sc <- argsDF$V2[which(argsDF$V1 == "fsc")]
if(is.null(argsL$hm)) print("pes projection file not supplied") else file.heatmap <- argsDF$V2[which(argsDF$V1 == "hm")]
if(is.null(argsL$o)) print("output not supplied") else output <- argsDF$V2[which(argsDF$V1 == "o")]

if(is.null(argsL$dpi) && is.null(argsL$ht) && is.null(argsL$wt)) { 
  print("no dimensions set: using standard values (height=7,width=11,dpi=300)") 
  ggheight <- 7
  ggwidth <- 11
  ggdpi <- 300
} else {
  ggheight <- argsDF$V2[which(argsDF$V1 == "ht")]
  ggwidth <- argsDF$V2[which(argsDF$V1 == "wt")]
  ggdpi <- argsDF$V2[which(argsDF$V1 == "dpi")]
}

#---------------------------------------------
#---------------- proper code ----------------

# function that extract the current directory and source doitRfunctions.R
locate.funtion.file <- function() {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
  function.file <- paste(sep="/", script.basename, "doitRfunctions.R")
  print(paste("Sourcing",function.file,"from",script.name))
  source(function.file)
}

# source all the functions & load the libraries
locate.funtion.file()
library("ggplot2"); library("gridExtra"); library("reshape2")

# plot simulation conditions
if (exists("file.potential") && exists("file.temperature") && exists("file.pressure") && exists("file.density") 
      && file.exists(paste(file.potential)) && file.exists(paste(file.temperature)) && file.exists(paste(file.pressure)) 
      && file.exists(paste(file.density))) { 
  editXVGfile(file.potential); editXVGfile(file.temperature); editXVGfile(file.pressure); editXVGfile(file.density)
  output.simcond <- paste0(output,"_simcond.png")
  gromacsSimCond(file.potential,file.temperature,file.pressure,file.density,output.simcond,ggheight,ggwidth,ggdpi)
}

# plot rmsd
if(exists("file.rmsd") && file.exists(paste(file.rmsd))) {
  editXVGfile(file.rmsd)
  ggxlab <- grep.labels(file.rmsd)[[1]]; ggylab <- grep.labels(file.rmsd)[[2]]
  ggtitle <- grep.title(file.rmsd)
  output.rmsd <- paste0(output,"_rmsd.png")
  generalplot(file.rmsd,output.rmsd,ggxlab,ggylab,ggtitle,ggheight,ggwidth,ggdpi)
}

# plot rgyr
if(exists("file.rgyr") && file.exists(paste(file.rgyr))) {
  editXVGfile(file.rgyr)
  ggxlab <- grep.labels(file.rgyr)[[1]]; ggylab <- grep.labels(file.rgyr)[[2]]
  ggtitle <- grep.title(file.rgyr)
  output.rgyr <- paste0(output,"_rgyr.png")
  generalplot(file.rgyr,output.rgyr,ggxlab,ggylab,ggtitle,ggheight,ggwidth,ggdpi)
}

# plot ss
if (exists("file.ss") && file.exists(paste(file.ss))) {
  editXVGfile(file.ss)
  ggxlab <- grep.labels(file.ss)[[1]]; ggylab <- grep.labels(file.ss)[[2]]
  ggtitle <- grep.title(file.ss)
  output.ss <- paste0(output,"_ss.png")
  gromacsSS(file.ss,output.ss,ggxlab,ggylab,ggtitle,ggheight,ggwidth,ggdpi)
}

# plot rmsf
if (exists("file.rmsf.bb") && exists("file.rmsf.sc") && file.exists(paste(file.rmsf.bb)) && file.exists(paste(file.rmsf.sc))) {
  editXVGfile(file.rmsf.bb); editXVGfile(file.rmsf.sc)
  ggxlab <- grep.labels(file.rmsf.bb)[[1]]; ggylab <- grep.labels(file.rmsf.bb)[[2]]
  ggtitle <- grep.title(file.rmsf.bb)
  output.rmsf <- paste0(output,"_rmsf.png")
  gromacsRMSF(file.rmsf.bb,file.rmsf.sc,output.rmsf,ggxlab,ggylab,ggtitle,ggheight,ggwidth,ggdpi)
}

# plot heat map
if (exists("file.heatmap") && file.exists(paste(file.heatmap))) {
  output.pes <- paste0(output,"_pes.png")
  gromacsPES(file.heatmap,output.pes,ggheight,ggwidth,ggdpi)
}

