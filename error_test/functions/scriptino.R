
args <- (commandArgs(T))

library(gdata);library(ggplot2);library(plotflow)
library(dplyr);library(grid)

binding.schem <- read.table(paste(args[1]),header=T,sep="",quote="",comment.char="#")
fluct.schem <- read.table(paste(args[2]),header=T,sep="",quote="",comment.char="#")
rmsf.data <- read.table(paste(args[3]),header=F,sep="",quote="")
names(rmsf.data) <- c("residue","nm")

bb.ggplot <- ggplot(rmsf.data,aes(x=rmsf.data$residue),environment = environment())

gg.rmsf <- bb.ggplot+ geom_line(aes(y=rmsf.data$nm)) + 
          labs(y = "nm", x = NULL, title = NULL) + 
          theme(legend.position="none",plot.margin=unit(c(0,1,0,1), "line"))

binding.sites <- ggplot(binding.schem,running=order_by(Domain, ~Start), aes(x=Type, ymin=Start, ymax=Stop, colour=Domain)) +
  scale_color_manual(values=c("white", "red","blue")) +
  geom_linerange(size = 5) + coord_flip() + theme_bw() + 
  labs(y = NULL, x = NULL, title = NULL) + theme(axis.line = element_blank(), 
  panel.grid.major = element_blank(),  
  panel.grid.minor = element_blank(), panel.border = element_blank(),   
  panel.background =  element_blank(), plot.background= element_blank(), axis.ticks = element_blank(),    
  legend.position ="none", axis.text = element_blank(), plot.margin=unit(c(0.5,1,0,1), "line"))

flucty.sites <- ggplot(fluct.schem, running=order_by(Domain, ~Start), aes(x=Type, ymin = Start, ymax = Stop, colour = Domain)) +
  scale_color_manual(values=c("white", "green3")) +
  geom_linerange(size = 5) + coord_flip() + theme_bw() + 
  labs(y = NULL, x = NULL, title = NULL) + theme(axis.line = element_blank(), 
  panel.grid.major = element_blank(),  
  panel.grid.minor = element_blank(), panel.border = element_blank(),   
  panel.background =  element_blank(), plot.background= element_blank(), axis.ticks = element_blank(),    
  legend.position ="none", axis.text = element_blank(), plot.margin=unit(c(-1.65,1,0,1), "line"))

gA <- ggplotGrob(gg.rmsf)
gB <- ggplotGrob(binding.sites)
gC <- ggplotGrob(flucty.sites)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
hello <- arrangeGrob(gA, gB,gC, ncol=1,top=textGrob("rmsf fluctuations", gp=gpar(fontsize=10),just="top"),heights = c(0.9,0.05,0.05))

ggsave(paste0(args[4],".png"),hello, height=7,width=11,dpi=300)
