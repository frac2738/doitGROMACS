
args <- (commandArgs(T))

# load libraries
library("ggplot2")
library("grid")
library("gridExtra")

# read the xvg file for the RMSD
hb.matrix <- read.table(paste(args[1]),header=F,sep=",",quote="")

hb.plot <-ggplot(hb.matrix, aes(x=V2,y=V5))
hb.plot + geom_point(aes(shape=V8,colour=V8,name="porcodio"))               +
  scale_shape(solid = FALSE)                                +
  coord_fixed(ratio = 1)+ggtitle("Hydrogen Bonds network")  + 
  ylab("residues")+xlab("residues")                 +
  theme(plot.title=element_text(size=10, vjust=2))  +
  theme(axis.title.x=element_text(size=9))          +
  theme(axis.title.y=element_text(size=9))          +
  theme(legend.title = element_blank())

ggsave(paste0(paste(args[1]),".png"),height=8,width=8, dpi=300)

