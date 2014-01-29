#Code to find missing jobs
#NWK 1-29-14
#Usage Rscript FindMissing.R Directory paramfile

library(plyr)
arguments <- commandArgs(trailingOnly=T)
allstats <- ldply(dir(arguments[1],pattern="stat",full.names=T),read.table,sep="\t",stringsAsFactors=F,header=F)

paramfile <- scan(arguments[2],what="character")
paramfile <- paramfile[grep("chunks",paramfile)]
nchunks <- Reduce("*",as.numeric(sapply(strsplit(paramfile,"="),"[",2)))-1

print(setdiff(0:nchunks,allstats[,1]))

