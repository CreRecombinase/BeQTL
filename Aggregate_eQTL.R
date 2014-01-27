#NWK 1/23/14
#Code to create Robj for eQTL results
#Usage: Rscript Aggregate_eQTL.R eqtl_directory eqtl_resultfile
#
library(plyr)
arguments <- commandArgs(trailingOnly=T)

directory <- arguments[1]
resultfile <- arguments[2]

eqtlfiles <- dir(directory,pattern="_eqtl")

alldata <- ldply(eqtlfiles,read.table,sep="\t",header=F,stringsAsFactors=F)

alldata <- alldata[!alldata$V1=="SNP",]
alldata <- unique(alldata)

colnames(alldata) <- c("SNP","Gene","t.Stat.Low","t.Stat.Med","t.Stat.High","Distance")

saveRDS(alldata,resultfile)
