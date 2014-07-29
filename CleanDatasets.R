library(sqldf)
library(jsonlite)


t.final <- function(x,n,tests,adjp){
  abs(p.adjust(pt(q = abs(x),df = n-2,lower.tail = F),method = "fdr",n = tests)-adjp)
}



#usage CleanDatasets.R type oldir annodir newdir snpcut genecut paramfile JobsizeInGB
args <- list()

data <- fromJSON(commandArgs(trailingOnly=T)[[1]])

paste0("Reading in snp file ",args$SNPfile)
snpdata <- read.csv.sql(data$inputsnpfile,sep="\t",header=T,eol="\n")
paste0("Reading in gene file ",args$GENEfile)
genedata <- read.csv.sql(data$inputexpfile,sep="\t",header=T,eol="\n")

paste0("Reading in SNP anno ",args$SNPanno)
paste0("Reading in gene anno ",args$Geneanno)


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

snpdata <- snpdata[snpcount>args$snpcut,]

genecount <- apply(genedata,1,function(x)sum(x>0))
genecount <- genecount/ncol(genedata);

genedata <- genedata[genecount>args$genecut,]

bsi <- ceiling(ncol(genedata)*log10(ncol(genedata)))
snpgenesize <- ceiling(sqrt((1024^3*args$GB)/(bsi*8)))

snpgenesize <- floor(snpgenesize/64)*64
snptotal <- nrow(snpdata)
snpchunks <- ceiling(snptotal/snpgenesize)
genetotal <- nrow(genedata)
genechunks <- ceiling(genetotal/snpgenesize)

t_thresh <- optimize(t.final,c(0,50),tol=0.00001,df=ncol(genedata)-2,tests=genetotal,adjp=0.05))

params <- c(snpfile=args$NewSNPfile,
            genefile=args$NewGENEfile,
            genecut=args$genecut,
            snpcut=args$snpcut,
            h5file=args$h5file,
            annofile=args$annofile,
            progfile=args$statfile,
            eqtlfile=args$eqtlfile,
            snpchunks=snpchunks,
            genechunks=genechunks,
            snptotal=snptotal,
            genetotal=genetotal,
            snpsize=snpgenesize,
            genesize=snpgenesize,
            casetotal=ncol(genedata),
            bsi=bsi,
            t_thresh=t_thresh,
            Streaming=args$Streaming,
            quantilenum=length(args$quantilepoints),
            quantilepoints=args$quantilepoints,
            cisdist=args$cisdist)

aparams <- paste(names(params),params,sep="=")

write(aparams,file=args$paramfile,sep="\n")


                  

write.table(snpdata,args$NewSNPfile,sep="\t",col.names=T,row.names=T,quote=F)
write.table(genedata,args$NewGENEfile,sep="\t",col.names=T,row.names=T,quote=F)


