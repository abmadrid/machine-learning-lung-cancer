# merge all files into a single data frame

library(plyr)
setwd(paste(getwd(),"/with_genes",sep=""))
file_list <- list.files()

mergeFiles<-function(file_list){
  dataset <- ldply(file_list, read.csv,header=FALSE)
  names(dataset)<-c("chromosome", "start", "end","num_probes","seg_mean","num_genes","most_freq" ,"cancer")
  #length(which(dataset$cancer==1))
  #length(which(dataset$cancer==-1))

  # write the data frame to file
  outfile<-"alldata_genes.csv"
  write.table(dataset,file=outfile,sep=",",row.names=FALSE)
}

mergeFiles(file_list)
