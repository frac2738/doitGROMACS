
#---------------- functions definition ----------------

#---------------- f() CheckExistance ----------------eckExistance <- function(a) {
checkExistance <- function(a) {
if (file.exists(a)) {
  cat('"',paste(a),'" -> ok\n')
  } else {
    cat('the file "',paste(a),'" does not exist\n')
  }
}

#---------------- f() edit xvg ----------------
# function that replaces all @ with # given a file
editXVGfile <- function(input) {
  for(f in input) {
    y <- gsub("@","#",readLines(f))
    cat(y,file=f,sep="\n") 
  }
}

#---------------- f() read.file ----------------
read.file <- function(filename) {
  input.table <- read.table(paste(filename),header=F,comment.char="#", sep="",quote="")
}

#---------------- f() grep.ss.types ----------------
grep.ss.types <- function(a) {
  lst <- Filter(function(u) grepl('^# s[0-9]+', u),readLines(paste(a)))
  return(gsub('.*\"(.*)\".*','\\1',lst))
}

#---------------- f() grep.labels ----------------
grep.labels <- function(a) {
  grep.xlst <- Filter(function(u) grepl('^#.*xaxis.*', u),readLines(paste(a)))
  grep.x <- gsub('.*\"(.*)\".*','\\1',grep.xlst)
  grep.ylst <- Filter(function(u) grepl('^#.*yaxis.*', u),readLines(paste(a)))
  grep.y <- gsub('.*\"(.*)\".*','\\1',grep.ylst)
  return(list(grep.x,grep.y))
}

#---------------- f() grep.title ----------------
grep.title <- function(a) {
  lst <- Filter(function(u) grepl('^#[ ]+title', u),readLines(paste(a)))
  return(gsub('.*\"(.*)\".*','\\1',lst))
}

#---------------- f() cbind.all ----------------
# function that adds "NA" to the shorter columns
cbind.all <- function (...) {
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(,n-nrow(x), ncol(x)))))
}

#---------------- f() generalplot(rmsd-rgyr) ----------------
generalplot <- function(a,o,x,y,t,h,w,dpi) {
  data.table <- read.file(a)
  names(data.table) <- c("time","value")
  data.ggplot <- ggplot(data.table,aes(x=time),environment=environment()) +
    xlab(x) + ylab(y) + ggtitle(t) + geom_line(aes(y=value)) +
    theme(legend.title=element_blank(), title = element_text(size=10))
  ggsave(o,height=h,width=w,dpi=dpi)
}

#---------------- f() generalComparisonPlot ----------------
generalComparisonPlot <- function(a,o,x,y,t) {
  data.ggplot <- ggplot(a,aes(x=time),environment=environment()) +
  xlab(x) + ylab(y) + ggtitle(t) + 
  geom_line(aes(y=value, colour=variable)) +
  theme(legend.title=element_blank(), title = element_text(size=10))
  ggsave(o,height=7,width=11,dpi=300) 
} 


# function used when there are multiple graphs on a plot
generalggplot <- function(a,x,y,t) {
  ggplot(a,aes(a$V1),environment = environment())+
  xlab(x) + ylab(y) + ggtitle(t) + geom_line(aes(y=a$V2))+
  theme(legend.title=element_blank(), title=element_text(size=10))
}

# ggplot if the data.set was melted with id=time
meltedggplot <- function(a,o,x,y,t,h,w,dpi) {
  ggplot(a,aes(time),environment = environment())+
  xlab(x) + ylab(y) + ggtitle(t) + geom_line(aes(y=value,colour=variable)) +
  theme(legend.title=element_blank(), title=element_text(size=10))
  ggsave(o,height=h,width=w,dpi=dpi)
}

gromacsSS <- function(a,o,x,y,t,h,w,dpi) {
  SStructure.table <- read.file(a)
  grep.names <- grep.ss.types(a)
  names(SStructure.table) <- c("time",grep.names)
  drops <- c("Chain_Separator")
  SStructure.table <- SStructure.table[,!names(SStructure.table) %in% drops]
  SStructure.table <- melt(SStructure.table,id.vars="time")
  ss.ggplot <- ggplot(SStructure.table,aes(x=SStructure.table$time),environment=environment())+
    xlab(x) + ylab(y) + ggtitle(t) +  
    geom_line(aes(y=SStructure.table$value,colour=SStructure.table$variable)) +
    theme(legend.title=element_blank(),title=element_text(size=10))
  ggsave(o,height=h,width=w,dpi=dpi)
}

#---------------- f() chain split ----------
# use flag with n° chain -ch=2 and then for (i in 1:$ch)
chain.split <- function(a) {
  rmsf.table <- data.frame()
  input.file <- read.file(a); names(input.file) <- c("residue","rmsf")
  chainA <- input.file[c(1:479),]; chainB <- input.file[c(480:958),]
  tmp.table <- cbind.all(rmsf.table, chainA); rmsf.table <- as.data.frame(tmp.table)
  tmp.table <- cbind.all(rmsf.table, chainB); rmsf.table <- as.data.frame(tmp.table)
  return(rmsf.table)
}

#---------------- f() RMSF ----------------
gromacsRMSF <- function(a,b,o,x,y,t,h,w,dpi) {
  rmsf.bb <- chain.split(a); rmsf.sc <- chain.split(b)
  rmsf.table <- cbind(rmsf.bb,rmsf.sc)
  names(rmsf.table) <- c("residue","bbA","residue2","bbB","residue3","scA","residue4","scB")

  bb.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment()) 
  bb.ggplot <- bb.ggplot + xlab(x) + ylab(y) +
    geom_line(aes(y=rmsf.table$bbA,colour="chain A")) +
    geom_line(aes(y=rmsf.table$bbB,colour="chain B")) + 
    theme(legend.position="none")
  
  sc.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment())
  sc.ggplot <- sc.ggplot + xlab(x) + ylab(y) +
    geom_line(aes(y=rmsf.table$scA,colour="chain A")) +
    geom_line(aes(y=rmsf.table$scB,colour="chain B")) + 
    theme(legend.position="none")
 
  savef <- arrangeGrob(bb.ggplot,sc.ggplot,ncol=1,top=textGrob(t, gp=gpar(fontsize=10),just="top"))
  ggsave(file=o,savef,height=h,width=w,dpi=dpi)  
}

gromacsRMSF2 <- function(a,b,o,x,y,t,h,w,dpi) {
  rmsf.bb <- chain.split(a); rmsf.sc <- chain.split(b)
  rmsf.table <- cbind(rmsf.bb,rmsf.sc)
  names(rmsf.table) <- c("residue","bbA","residue2","bbB","residue3","scA","residue4","scB")

  bb.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment()) 
  bb.ggplot <- bb.ggplot + xlab(x) + ylab(y) + theme(axis.title.x = element_blank()) +
    geom_line(aes(y=rmsf.table$bbA,colour="chain A")) + 
    geom_line(aes(y=rmsf.table$bbB,colour="chain B")) + 
    theme(legend.position="none")
  
  sc.ggplot <- ggplot(rmsf.table,aes(x=rmsf.table$residue),environment = environment())
  sc.ggplot <- sc.ggplot + xlab(x) + ylab(y) +
    geom_line(aes(y=rmsf.table$scA,colour="chain A")) +
    geom_line(aes(y=rmsf.table$scB,colour="chain B")) + 
    theme(legend.position="none")

  binding.schem <- read.table("/export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/schematics_bindingSites",header=T,sep="",quote="",comment.char="#")
  fluct.schem <- read.table("/export/francesco/gitHub/doitGROMACS/doitgromacs-v.1.6/schematics_motileRegions",header=T,sep="",quote="",comment.char="#")

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
    legend.position ="none", axis.text = element_blank(), plot.margin=unit(c(-1.6,1,0,1), "line"))

  gA  <- ggplotGrob(bb.ggplot)
  gAA <- ggplotGrob(sc.ggplot)
  gB  <- ggplotGrob(binding.sites)
  gC  <- ggplotGrob(flucty.sites)
  maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gAA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  gC$widths[2:5] <- as.list(maxWidth)
  savef <- arrangeGrob(gA,gB,gC,gAA, ncol=1,top=textGrob("rmsf fluctuations", gp=gpar(fontsize=10),just="top"),heights = c(0.45,0.05,0.05,0.45))

  #savef <- arrangeGrob(bb.ggplot,sc.ggplot,ncol=1,top=textGrob(t, gp=gpar(fontsize=10),just="top"))
  ggsave(file=o,savef,height=h,width=w,dpi=dpi)  
}

#---------------- f() coordinate extractor ----------------
coord.extractor <- function(l,i) {
  coord.vector <- l[i,]; names(coord.vector) <- NULL
  coord.vector <- unlist(c(coord.vector)); coord.vector <- coord.vector[1:length(coord.vector)-1]
  return(coord.vector)
}

#---------------- f() HEAT MAP ----------------
gromacsPES <- function(a,o,h,w,dpi) {
  heatdata <- read.table(paste(a),header=F,sep="",quote="",blank.lines.skip=F)
  x.coordinate <- coord.extractor(heatdata,1)
  y.coordinate <- coord.extractor(heatdata,2)
  
  energy.table <- heatdata[(3:34),] 
  names(energy.table) <- x.coordinate
  energy.table[is.na(energy.table)] <- max(energy.table[,1])
  energy.table$y <- y.coordinate
  longenergy <- melt(energy.table,id="y")
  
  ggheat <- ggplot(longenergy,aes(x=longenergy$variable,y=longenergy$y, fill=value),environment = environment())
  ggheat <- ggheat + geom_tile()+ theme(axis.text.x=element_blank(),axis.text.y=element_blank()) + 
    xlab("pc1") + ylab("pc2") + ggtitle("Gibbs Energy Landscape [kJ/mol]") 
    #scale_fill_gradient(low="blue",high="red")
  ggsave(o,height=h,width=w,dpi=dpi)
}

gromacsSimCond <- function(a,b,c,d,o,h,w,dpi) {
  potential <- read.file(a); temperature <- read.file(b)
  pressure <- read.file(c); density <- read.file(d)

  potxlab <- grep.labels(file.potential)[[1]]; potylab <- grep.labels(file.potential)[[2]]
  pottitle <- "potential"
  pot.ggplot <- generalggplot(potential,potxlab,potylab,pottitle)
  tempxlab <- grep.labels(file.temperature)[[1]]; tempylab <- grep.labels(file.temperature)[[2]]
  temptitle <- "temperature"
  temp.ggplot <- generalggplot(temperature,tempxlab,tempylab,temptitle)
  pressxlab <- grep.labels(file.pressure)[[1]]; pressylab <- grep.labels(file.pressure)[[2]]
  presstitle <- "pressure"
  press.ggplot <- generalggplot(pressure,pressxlab,pressylab,presstitle)
  denxlab <- grep.labels(file.density)[[1]]; denylab <- grep.labels(file.density)[[2]]
  dentitle <- "density"
  density.ggplot <- generalggplot(density,denxlab,denylab,dentitle)

  ggtitle <- grep.title(file.potential)
  savef <- arrangeGrob(pot.ggplot,temp.ggplot,press.ggplot,density.ggplot,ncol=2, 
    top=textGrob(ggtitle, gp=gpar(fontsize=10),just="top"))
  grid.draw(savef)
  ggsave(o,savef ,height=h,width=w,dpi=dpi)
}

gromacsSAS <- function(a,o,x,y,t,h,w,dpi) {
  data.set <- read.file(a)
  names(data.set) <- c("time","Hydrophobic","Hydrophilic","Total","D Gsolv")
  data.set.long <- melt(data.set,id="time")
  meltedggplot(data.set.long,o,x,y,t,h,w,dpi) 
}

