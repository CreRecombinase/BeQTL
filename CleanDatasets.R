library(sqldf)
#usage CleanDatasets.R SNPfile GENEfile SNPanno GENEanno NewSNPfile NewGENEfile
oargs <- commandArgs(trailingOnly=T)
args <- list()

args$SNPfile <- oargs[[1]]
args$GENEfile <- oargs[[2]]
args$SNPanno <- oargs[[3]]
args$Geneanno <- oargs[[4]]
args$NewSNPfile <- oargs[[5]]
args$NewGENEfile <- oargs[[6]]

snpdata <- read.csv.sql(args$SNPfile,sep="\t",header=T,eol="\n")
genedata <- read.csv.sql(args$GENEfile,sep="\t",header=T,eol="\n")

snpanno <- read.csv.sql(args$SNPanno,sep="\t",header=F,eol="\n")
geneanno <- read.csv.sql(args$Geneanno,sep="\t",header=F,eol="\n")


snpdata <- snpdata[!duplicated(snpdata[,1]),]
genedata <- genedata[!duplicated(genedata[,1]),]

rownames(snpdata) <- snpdata[,1]
rownames(genedata) <- genedata[,1]

snpdata <- snpdata[,-1]
genedata <- genedata[,-1]

colnames(snpdata) <- substr(colnames(snpdata),1,12)
colnames(genedata) <- substr(colnames(genedata),1,12)

snpdata <-  snpdata[,colnames(snpdata) %in% colnames(genedata)]
genedata <- genedata[,colnames(snpdata)]

snpanno <- snpanno[snpanno[,1] %in% rownames(snpdata),]
geneanno <- geneanno[geneanno[,1] %in% rownames(genedata),]

snpdata <- snpdata[snpanno[,1],]
geneanno <- genedata[geneanno[,1],]

write.table(snpdata,args$NewSNPfile,sep="\t",col.names=T,row.names=T,quote=F)
write.table(genedata,args$NewGENEfile,sep="\t",col.names=T,row.names=T,quote=F)


