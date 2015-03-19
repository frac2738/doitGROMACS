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
#         - all the xvg files have been generated using GROMACS
#
#         If a system differ for all, of some of the conditions listed above,
#         wait...I'll release a better version at some point in the future.
#
#------------------------------------------------------------------------------
#
# Usage:
# Rscript doitRGROMACS.R [-x=xx] [-t=xx] [-u=xx] [-p=xx] [-d=xx] [-g=xx]
#                  [-ss=xx] [-fb=xx] [-fsc=xx] [-help]
#
# -help   -  ask for help
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

      -help   -  ask for help
      -u      -  read potential energy file (.xvg)
      -t      -  read temperature file (.xvg)
      -p      -  read pressure file (.xvg)
      -x      -  read density file (.xvg)
      -d      -  read rmsd file (.xvg)
      -g      -  read gyration radius file (.xvg)
      -ss     -  read secondary structure file (.xvg)
      -fb     -  read rmsf file [backbone] (.xvg)
      -fsc    -  read rmsf file [side chains] (.xvg)
      -hm     -  read pes projection (heat map)

      NOTE1 :  flags x,t,u,p are read together, the absence of one (or more) of
               them will cause the function to exit without plotting anything. 
      
")
   q(save="no")
}

## Parse the arguments in the form of -arg=value
parseArgs <- function(x) strsplit (sub("^-","",x), "=")   # split flags from values
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
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


#---------------------------------------------
#---------------- Proper code ----------------

# load libraries
library("ggplot2")
library("gridExtra")
library("reshape2")

#---------------- f() edit xvg ----------------
# function that replaces all @ with # given a file
editXVGfile <- function(input) {
  for(f in input) {
    x <- readLines(f)
    y <- gsub("@","#",x)
    cat(y,file=f,sep="\n") 
  }
}

#---------------- f() cbind.all ----------------
cbind.all <- function (...) 
{
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - 
    nrow(x), ncol(x)))))
}

#---------------- f() RMSD ----------------
gromacsRMSD <- function(a) {
  rmsd.table <- read.table(paste(a),header=F,comment.char="#", sep="",quote="")
  names(rmsd.table) <- c("time","rmsd")
  rmsd.ggplot <- ggplot(rmsd.table,aes(x=rmsd.table$time),environment=environment()) +
    xlab("[ns]") + ylab("[nm]") + ggtitle("rmsd") + geom_line(aes(y=rmsd.table$rmsd))  +
    theme(legend.title=element_blank(), title = element_text(size = rel(2)))
  ggsave("rmsd.png",height=7,width=11,dpi=300)
}

if(exists("file.rmsd") && file.exists(paste(file.rmsd))) {
  editXVGfile(file.rmsd)
  gromacsRMSD(file.rmsd)
}

#---------------- f() Gyration radius ----------------
gromacsRGYRATION <- function(a) {
  rgyr.table <- read.table(paste(a),header=F,comment.char="#",sep="",quote="")
  names(rgyr.table) <- c("time","rgyr")
  rgyr.ggplot <- ggplot(rgyr.table,aes(x=rgyr.table$time),environment=environment()) +
    xlab("[ps]") + ylab("[nm]") + ggtitle("gyration radius") + geom_line(aes(y=rgyr.table$rgyr)) + 
    theme(legend.title=element_blank(), title = element_text(size = rel(2)))
  ggsave("rgyr.png",height=7,width=11,dpi=300)
}
 
if(exists("file.rgyr") && file.exists(paste(file.rgyr))) {
  editXVGfile(file.rgyr)
  gromacsRGYRATION(file.rgyr)
}

#---------------- f() Secondary Structure ----------------
gromacsSS <- function(a) {
  SStructure.table <- read.table(paste(a),header=F,comment.char="#",sep="",quote="")
  list.names <- list()
  for(i in 1:ncol(SStructure.table)) {
    list.names[i] <- i
  }
  # try also con grep to grep the names

  #names(SStructure.table) <- list.names
  names(SStructure.table) <- c("time","Structure","Coil","B-Sheet","B-Bridge",
    "Bend","Turn","A-Helix","5-Helix","3-Helix","Chain_Separator")
  drops <- c("Chain_Separator")
  SStructure.table <- SStructure.table[,!names(SStructure.table) %in% drops]
  SStructure.table <- melt(SStructure.table,id.vars="time")

  ss.ggplot <- ggplot(SStructure.table,aes(x=SStructure.table$time),environment=environment())+
    xlab("[ns]") + ylab("# residues") + ggtitle("SS %") +  
    geom_line(aes(y=SStructure.table$value,colour=SStructure.table$variable)) +
    theme(legend.title=element_blank())
  ggsave("ss.png",height=7,width=11,dpi=300)
}

if (exists("file.ss") && file.exists(paste(file.ss))) {
   editXVGfile(file.ss)
   gromacsSS(file.ss)
}

#---------------- f() Simulation conditions ----------------

gromacsSimCond <- function(a,b,c,d) {
  potential <- read.table(paste(a), comment.char="#",header=F,sep="",quote="")
  temperature <- read.table(paste(b), comment.char="#",header=F,sep="",quote="")
  pressure <- read.table(paste(c), comment.char="#",header=F,sep="",quote="")
  density <- read.table(paste(d), comment.char="#",header=F,sep="",quote="")

  Simtime <- potential$V1/1000 
  pot.ggplot <- ggplot(potential, aes(x=Simtime),environment = environment()) +
    xlab("[ns]") + ylab("[kJ/mol]") + ggtitle("Potential") + geom_line(aes(y=potential$V2)) + 
    theme(axis.title.x=element_blank(),legend.title=element_blank(), title=element_text(size = rel(2)))

  temp.ggplot <- ggplot(temperature, aes(x=Simtime),environment = environment()) +
    xlab("[ns]") + ylab("[K]") + ggtitle("Temperature") + geom_line(aes(y=temperature$V2)) + 
    theme(axis.title.x=element_blank(),legend.title=element_blank(), title=element_text(size = rel(2)))
  press.ggplot <- ggplot(pressure, aes(x=Simtime),environment = environment()) +
    xlab("[ns]") + ylab("[bar]") + ggtitle("Pressure") + geom_line(aes(y=pressure$V2)) + 
    theme(legend.title=element_blank(), title=element_text(size = rel(2)))
  desity.ggplot <- ggplot(density, aes(x=Simtime),environment = environment()) +
    xlab("[ns]") + ylab("[kg/m^3]") + ggtitle("density") +
    geom_line(aes(y=density$V2)) + theme(legend.title=element_blank(), title=element_text(size = rel(2)))

  savef <- arrangeGrob(pot.ggplot,temp.ggplot,press.ggplot,desity.ggplot,ncol=2, main=textGrob("Simulation Conditions", gp=gpar(fontsize=85)))
  ggsave(file="simcond.png",savef ,height=7,width=12,dpi=300)
}

if (exists("file.potential") && exists("file.temperature") && exists("file.pressure") && exists("file.density") 
      && file.exists("file.potential") && file.exists("file.temperature") && file.exists("file.pressure") 
      && exists("file.density")) { 
  editXVGfile(file.potential)
  editXVGfile(file.temperature)
  editXVGfile(file.pressure)
  editXVGfile(file.density)
  gromacsSimCond(file.potential,file.temperature,file.pressure,file.density)
}

#---------------- f() RMSF ----------------
gromacsRMSF <- function(a,b) {
  rmsf.table <- data.frame()
 
  input.file <- read.table(paste(a),comment.char="#",header=F,sep="",quote="")
  names(input.file) <- c("residue","rmsf")
  chainA <- input.file[c(1:479),]
  chainB <- input.file[c(480:958),]
  tmp.table <- cbind.all(rmsf.table, chainA)
  rmsf.table <- as.data.frame(tmp.table)
  tmp.table <- cbind.all(rmsf.table, chainB)
  rmsf.table <- as.data.frame(tmp.table)
  
  input2.file <- read.table(paste(b),comment.char="#",header=F,sep="",quote="")
  names(input2.file) <- c("residue","rmsf")
  chainA <- input2.file[c(1:479),]
  chainB <- input2.file[c(480:958),]
  tmp2.table <- cbind.all(rmsf.table, chainA)
  rmsf.table <- as.data.frame(tmp2.table)
  tmp2.table <- cbind.all(rmsf.table, chainB)
  rmsf.table <- as.data.frame(tmp2.table)

  names(rmsf.table) <- c("residue","bbA","residue2","bbB","residue3","scA","residue4","scB")

  bb.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment()) 
  bb.ggplot <- bb.ggplot + xlab("residue") + ylab("[nm]") + ggtitle("Backbone") +
    geom_line(aes(y=rmsf.table$bbA,colour="chain A")) +
    geom_line(aes(y=rmsf.table$bbB,colour="chain B")) + 
    theme(legend.position="none")

  sc.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment())
  sc.ggplot <- sc.ggplot + xlab("residue") + ylab("[nm]") + ggtitle("Side Chain") +
    geom_line(aes(y=rmsf.table$scA,colour="chain A")) +
    geom_line(aes(y=rmsf.table$scB,colour="chain B")) + 
    theme(legend.position="none")

 
  savef <- arrangeGrob(bb.ggplot,sc.ggplot,ncol=1,main=textGrob("RMSF", gp=gpar(fontsize=85)) )
  ggsave(file="rmsf.png",savef,height=7,width=11,dpi=300)   
}

if (exists("file.rmsf.bb") && exists("file.rmsf.sc") && file.exists(paste(file.rmsf.bb)) && file.exists(paste(file.rmsf.sc))) {
  editXVGfile(file.rmsf.bb)
  editXVGfile(file.rmsf.sc)
  gromacsRMSF(file.rmsf.bb,file.rmsf.sc)
}

#---------------- f() HEAT MAP ----------------

gromacsPES <- function(a) {
  heatdata <- read.table(paste(a),header=F,sep="",quote="",blank.lines.skip=F)
  x.coordinate <- heatdata[1,]
  names(x.coordinate) <- NULL
  x.coordinate <- unlist(c(x.coordinate))
  x.coordinate <- x.coordinate[1:length(x.coordinate)-1]
  
  y.coordinate <- heatdata[2,]
  names(y.coordinate) <- NULL
  y.coordinate <- unlist(c(y.coordinate))
  y.coordinate <- y.coordinate[1:length(y.coordinate)-1]

  energy.table <- heatdata[(3:34),] 
  names(energy.table) <- x.coordinate
  energy.table[is.na(energy.table)] <- max(energy.table[,1])
  energy.table$y <- y.coordinate
  longenergy <- melt(energy.table,id="y")
  
  ggheat <- ggplot(longenergy,aes(x=longenergy$variable,y=longenergy$y, fill=value),environment = environment())
  ggheat <- ggheat + geom_tile()+ theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + 
    xlab("pc1") + ylab("pc2") + ggtitle("Gibbs Energy Landscape [kJ/mol]") 
    #scale_fill_gradient(low="blue",high="red")
  ggsave(file="pes.png",height=7,width=11,dpi=300)
}

if (exists("file.heatmap") && file.exists(paste(file.heatmap))) {
  gromacsPES(file.heatmap)
}

