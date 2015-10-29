# load libraries
# Rscipt script.R output file1 file2 file3 ...

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

args <- (commandArgs(TRUE))

data.table <- data.frame()          # initialise an empty data.frame

for (i in 2:length(args)) {
  input.file <- read.file(args[i])
  names(input.file) <- c("time",gsub('.xvg',"\\1",args[i]))
  tmp.table <- cbind.all(data.table, input.file)
  data.table <- as.data.frame(tmp.table)
}

ggxlab <- grep.labels(args[2])[[1]]; ggylab <- grep.labels(args[2])[[2]]
ggtitle <- grep.title(args[2])
output <- args[1]
output.name <- paste0(output,".png",sep="")

data.table <- melt(data.table, id="time")

generalComparisonPlot(data.table,output.name,ggxlab,ggylab,ggtitle)

