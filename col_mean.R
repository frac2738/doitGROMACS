args <- (commandArgs(TRUE))

parseArgs <- function(x) strsplit (sub("^-","",x), "=")   # split flags from values
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- (argsDF$V1)

if(!is.null(argsL$sas)) file.values <- argsDF$V2[which(argsDF$V1 == "sas")]
if(!is.null(argsL$hb)) file.values <- argsDF$V2[which(argsDF$V1 == "hb")]

values <- read.table(paste(file.values),header=F,comment.char="#",sep="",quote="")

## default behaviours
if(!is.null(argsL$sas)) {
  sas_mean <- mean(values$V4)
  print(sas_mean)
}

if(!is.null(argsL$hb)) {
  hb_mean <- mean(values$V2)
  print(hb_mean)
}
