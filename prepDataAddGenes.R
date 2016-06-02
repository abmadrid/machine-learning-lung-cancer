# data preprocessing of files:
#   convert absolute copy number
#   add number of genes, most frequent genes, cancer label
setwd(paste(getwd(),"/raw_data",sep=""))
raw.files<-list.files()
raw.files<-paste(getwd(),"/",raw.files,sep="")

setwd("../gene_data")
gene.files<-list.files()
gene.files<-paste(getwd(),"/",gene.files,sep="")

all_files<-cbind(raw.files,gene.files)

addGenes<-function(all_files){
  for (i in 1:dim(all_files)[1]){
    print(paste("processing file # ",i))

    infile_raw<-all_files[i,1]
    infile_genes<-all_files[i,2]
    data.raw<-read.table(infile_raw,skip=1)
    data<-data.raw[,2:6]
    names(data)<-c("chromosome", "start", "end", "num_probes", "seg_mean")
    data.genes<-read.table(infile_genes,skip=1)
    names(data.genes)<-c("chromosome","start","end","gene")

    # absolute copy number conversion (2^seg_mean)*2
    data$seg_mean<-(2^data$seg_mean)*2
    data$chr<-data$chromosome
    data$num_genes<-rep(NA,dim(data)[1])
    data$most_freq<-rep(NA,dim(data)[1])

    # add genes
    for (j in 1:dim(data)[1]){
      c<-paste("chr",data$chr[j],sep="")
      s<-data$start[j]
      e<-data$end[j]
      data.g.chr<-data.genes[which(data.genes$chromosome==c),]
      data$num_genes[j]<-length(which((data.g.chr$start>=s)&(data.g.chr$end<=e)))
      g<-as.character(data.genes$gene[which((data.g.chr$start>=s)&(data.g.chr$end<=e))])
      data$most_freq[j]<-sort(table(g),decreasing = TRUE)[1]
    }
    data$most_freq[is.na(data$most_freq)]<-0

    # append labels for cancer/no cancer
    if (grepl("nocnv",infile_raw)) { 
      data$cancer<-rep(-1,dim(data)[1]) 
    } 
    else { data$cancer<-rep(1,dim(data)[1]) }

    # save to outfile
    outfile<-paste("../with_genes/",i,".csv",sep="")
    write.table(data[,-6],file=outfile,sep=",",row.names=FALSE,col.names=FALSE)
  }

  print("done")
}

addGenes(all_files)
