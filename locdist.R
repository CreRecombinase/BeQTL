#Script for finding distant and local eQTL
#Usage Rscript locdist.R snpanno.rds geneanno.rds eqtl.rds 
arguments <- commandArgs(trailingOnly=T)

snpAnno <- readRDS(arguments[1])
geneAnno <- readRDS(arguments[2])

eqtlData <- readRDS(arguments[3])

eqtlData$Gene <- gsub("(.+)\\|.+","\\1",eqtlData$Gene)

allSnps <- snpAnno[eqtlData$SNP,]
allGenes <- geneAnno[eqtlData$Gene,]

sameChromosome <- allSnps$Chrom==allGenes$Chrom

sameChromosome[is.na(sameChromosome)] <- FALSE

cisSnps <- allSnps[sameChromosome,]
cisGenes <- allGenes[sameChromosome,]

cisDist <- pmin(abs(cisSnps$Position-cisGenes$Start),abs(cisSnps$Position-cisGenes$Stop))
eqtlData[,"Distance"]<- as.numeric(eqtlData$Distance)
eqtlData[sameChromosome,"Distance"]<-cisDist
saveRDS(eqtlData,arguments[3])
