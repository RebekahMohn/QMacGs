
for(i in c(1:12)){
  
  position<-read.delim(paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/MAC99_MACREF_HQSNP_C",i,"_maf01.map",sep=""),header=FALSE)
  position$V2<-paste("V",row.names(position),sep="")
  write.table(position,paste("C:/Users/rmohn/Desktop/10_Analysis/114_GEA/lfmmOut/MAC99_MACREF_HQSNP_C",i,"_maf01.map2",sep=""),sep="\t",row.names = F,quote = F, col.names = F)
  
  }