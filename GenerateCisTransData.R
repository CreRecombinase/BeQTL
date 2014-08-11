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
  Gene_anno[allGenedata[[i]],"Cancer"] <- bitwOr(snpnumber[[i]],Gene_anno[allGenedata[[i]],"Cancer"])
}

SNPcan <- matrix(0,nrow=nrow(SNP_anno),ncol=length(snpnumber))
Genecan <- matrix(0,nrow=nrow(Gene_anno),ncol=length(snpnumber))
rownames(SNPcan) <- rownames(SNP_anno)
rownames(Genecan) <- rownames(Gene_anno)
colnames(SNPcan) <- names(snpnumber)
colnames(Genecan) <- names(snpnumber)

SNPcl <- split(SNP_anno[SNP_anno$Cancer!=0,],f = SNP_anno[SNP_anno$Cancer!=0,"Chromosome"])
Genecl <- split(Gene_anno[Gene_anno$Cancer!=0,],f=Gene_anno[Gene_anno$Cancer!=0,"Chromosome"])

for(i in 1:length(SNPcl)){
  print(i)
  startmat <- outer(SNPcl[[i]][,"Position"],Genecl[[i]][,"StartPosition"],function(x,y)(abs(x-y)<1e6))
  stopmat <- outer(SNPcl[[i]][,"Position"],Genecl[[i]][,"StopPosition"],function(x,y)(abs(x-y)<1e6))
  CancerMat <- outer(SNPcl[[i]][,"Cancer"],Genecl[[i]][,"Cancer"],bitwAnd)
  CancerMat[!(startmat|stopmat)]<-0
  
  
  CancCount <- sapply(snpnumber,function(x)matrix(bitwAnd(x,CancerMat)/x,nrow = nrow(CancerMat),ncol=ncol(CancerMat)))
  for(j in 1:length(snpnumber)){
    paste0("Inner Loop:" ,j)
    SNPcan[SNPcl[[i]][,"RS.ID"],j] <- apply(CancerMat,1,function(x)sum(bitwAnd(snpnumber[j],x)/snpnumber[j]))
    Genecan[Genecl[[i]][,"GeneName"],j] <- apply(CancerMat,2,function(x)sum(bitwAnd(snpnumber[j],x)/snpnumber[j]))
  }
}


SNP_anno <- data.frame(SNP_anno,SNPcan,stringsAsFactors=F)
Gene_anno <- data.frame(Gene_anno,Genecan,stringsAsFactors=F)
saveRDS(SNP_anno,"/groups/becklab/BeQTL/SNP_anno_ct.rds")
write.table(SNP_anno,"/groups/becklab/BeQTL/SNP_anno_ct.txt",sep="\t",col.names=T,row.names=F)
saveRDS(Gene_anno,"/groups/becklab/BeQTL/Gene_anno_ct.rds")
write.table(Gene_anno,"/groups/becklab/BeQTL/Gene_anno_ct.rds",sep="\t",col.names=T,row.names=F)

