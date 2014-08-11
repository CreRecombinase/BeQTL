#Code to generate CancerID assignment for each SNP and Gene

readData <- function(x){
  data <- system(paste0("cut -f 1 ",x),intern=T)[-1]
}
allSNPfiles <- dir("/groups/becklab/BeQTL/Inputs",pattern="snp_[A-Za-z]+.txt",full.names=T)
allSNPdata <- lapply(allSNPfiles,FUN=readData)
names(allSNPdata) <- tolower(gsub(".+snp_(.+).txt","\\1",allSNPfiles))
snpnumber <- 2^(0:(length(allSNPfiles)-1))
names(snpnumber) <- names(allSNPdata)
SNP_anno <- read.table("~/snpanno.txt",sep="\t",header=T,stringsAsFactors = F)
rownames(SNP_anno) <- SNP_anno$RS.ID
SNP_anno$Cancer <- 0

for(i in names(snpnumber)){
  print(i)
  SNP_anno[allSNPdata[[i]],"Cancer"] <- bitwOr(snpnumber[[i]],SNP_anno[allSNPdata[[i]],"Cancer"])
}

allGenefiles <- dir("/groups/becklab/BeQTL/Inputs",pattern="seq_",full.names=T)
allGenedata <- lapply(allGenefiles,FUN=readData)
names(allGenedata) <- names(allSNPdata)
Gene_anno <- read.table("~/geneanno.txt",sep="\t",header=F,stringsAsFactors=F)
colnames(Gene_anno) <- c("GeneName","Chromosome","StartPosition","StopPosition","GID")
rownames(Gene_anno) <- Gene_anno$GeneName
Gene_anno$Cancer <- 0

for(i in names(snpnumber)){
  print(i)
  Gene_anno[allGenedata[[i]],"Cancer"] <- bitwOr(snpnumber[[i]],Gene_anno[allGenedata[[i]],"Cancer"])
}

SNP_anno$CisGenes <- 0
Gene_anno$CisSNPs <- 0
SNPcl <- split(SNP_anno[SNP_anno$Cancer!=0,],f = SNP_anno[SNP_anno$Cancer!=0,"Chromosome"])
Genecl <- split(Gene_anno[Gene_anno$Cancer!=0,],f=Gene_anno[Gene_anno$Cancer!=0,"Chromosome"])

for(i in 1:length(SNPcl)){
  print(i)
  startmat <- outer(SNPcl[[i]][,"Position"],Genecl[[i]][,"StartPosition"],function(x,y)(abs(x-y)<1e6))
  stopmat <- outer(SNPcl[[i]][,"Position"],Genecl[[i]][,"StopPosition"],function(x,y)(abs(x-y)<1e6))
  starmat <- startmat|stopmat
  SNP_anno[SNPcl[[i]][,"RS.ID"],"CisGenes"] <- rowSums(startmat)
  Gene_anno[Genecl[[i]][,"GeneName"],"CisSNPs"] <- colSums(startmat))
}

saveRDS(SNP_anno,"/groups/becklab/BeQTL/SNP_anno_ct.rds")
write.table(SNP_anno,"/groups/becklab/BeQTL/SNP_anno_ct.txt",sep="\t",col.names=T,row.names=F)
saveRDS(Gene_anno,"/groups/becklab/BeQTL/Gene_anno_ct.rds")
write.table(Gene_anno,"/groups/becklab/BeQTL/Gene_anno_ct.rds",sep="\t",col.names=T,row.names=F)

