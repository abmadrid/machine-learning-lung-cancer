# parsing original data files for just the chromosome #, start position, and end position
# from working dir
setwd(paste(getwd(),"/raw_data",sep=""))
files<-list.files()

parsePos<-function(files){
  for (i in 1:length(files)){
    infile<-files[i]
    data.all<-read.table(infile,skip=1)
    names(data.all)<-c("sample", "chromosome", "start", "end", "num_probes", "seg_mean")
    # subset of data that contains chromosome number, start position, end position
    data.pos<-data.all[2:4]
    data.pos$chromosome<-paste("chr",data.pos$chromosome,sep="")
    # write chromosome number, start position, end position to outfile in uscs dir
    setwd("../for_uscs")
    outfile<-paste(as.character(i),"_pos.txt",sep="")
    write.table(data.pos,file=outfile,sep=" ",row.names=FALSE,quote=FALSE,col.names = FALSE)
    setwd("../raw_data")
  }
  print("done")
}
parsePos(files)
