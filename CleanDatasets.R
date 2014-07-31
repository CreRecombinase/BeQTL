library(sqldf)
library(jsonlite)


t.final <- function(x,n,tests,adjp){
  abs(p.adjust(pt(q = abs(x),df = n-2,lower.tail = F),method = "fdr",n = tests)-adjp)
}



#usage CleanDatasets.R type oldir annodir newdir snpcut genecut paramfile JobsizeInGB
args <- list()

data <- fromJSON(commandArgs(trailingOnly=T)[[1]])

paste0("Reading in snp file ",data$inputsnpfile)
snpdata <- read.csv.sql(data$inputsnpfile,sep="\t",header=T,eol="\n")
paste0("Reading in gene file ",data$inputgenefile)
genedata <- read.csv.sql(data$inputgenefile,sep="\t",header=T,eol="\n")


snpdata <- snpdata[!duplicated(snpdata[,1]),]
genedata <- genedata[!duplicated(genedata[,1]),]

rownames(snpdata) <- snpdata[,1]
rownames(genedata) <- genedata[,1]

snpdata <- snpdata[,-1]
genedata <- genedata[,-1]


colnames(snpdata) <- substr(colnames(snpdata),1,12)
paste0("first snpcols: ",Reduce(paste,head(colnames(snpdata))))
colnames(genedata) <- substr(colnames(genedata),1,12)
paste0("first genecols: ",Reduce(paste,head(colnames(genedata))))

paste0("Subsetting snpdata")
snpdata <-  snpdata[,colnames(snpdata) %in% colnames(genedata)]
paste0("subsetting genedata")
genedata <- genedata[,colnames(snpdata)]


snpcount <- apply(snpdata,1,function(x)sum(sort(tabulate(x+1),decreasing=T)[-1]))
snpcount <- snpcount/ncol(snpdata)

snpdata <- snpdata[snpcount>data$snpcut,]

genecount <- apply(genedata,1,function(x)sum(x>0))
genecount <- genecount/ncol(genedata);

genedata <- genedata[genecount>data$genecut,]

bsi <- ceiling(ncol(genedata)*log10(ncol(genedata)))
if(data$Streaming){
  snpgenesize <- ceiling(sqrt((1024^3*data$GB)/8))
}else{
  snpgenesize <- ceiling(sqrt((1024^3*data$GB)/(bsi*8)))
}
snpgenesize <- floor(snpgenesize/64)*64
snptotal <- nrow(snpdata)
snpchunks <- ceiling(snptotal/snpgenesize)
genetotal <- nrow(genedata)
genechunks <- ceiling(genetotal/snpgenesize)

t_thresh <- optimize(t.final,c(0,10),tol=0.00001,n=ncol(genedata),tests=genetotal,adjp=0.05,maximum=F)[["minimum"]]

params <- c(snpfile=data$outputsnpfile,
            genefile=data$outputgenefile,
            genecut=data$genecut,
            snpcut=data$snpcut,
            h5file=data$h5file,
            annofile=data$annofile,
            progfile=data$statfile,
            eqtlfile=data$eqtlfile,
            snpchunks=snpchunks,
            genechunks=genechunks,
            snptotal=snptotal,
            genetotal=genetotal,
            snpsize=snpgenesize,
            genesize=snpgenesize,
            casetotal=ncol(genedata),
            bsi=bsi,
            t_thresh=t_thresh,
            Streaming=data$Streaming,
            quantilenum=length(data$quantilepoints),
            quantilepoints=data$quantilepoints,
            cisdist=data$cisdist)

aparams <- paste(names(params),params,sep="=")

write(aparams,file=data$outputfile,sep="\n")


                  

write.table(snpdata,args$NewSNPfile,sep="\t",col.names=T,row.names=T,quote=F)
write.table(genedata,args$NewGENEfile,sep="\t",col.names=T,row.names=T,quote=F)


