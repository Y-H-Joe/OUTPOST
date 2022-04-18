args=commandArgs(T)
counts_dp=args[1]
kaiju_dp=args[2]
output=args[3]

counts=read.csv(counts_dp,header = FALSE,sep="\t")
kaiju=read.csv(kaiju_dp,header=FALSE,sep = "\t")

sample_num=dim(counts)[2]-1
withTaxaID=TRUE


kaiju_counts=merge(kaiju,counts,by="V1")
if (withTaxaID==TRUE){
  kaiju_counts$sum <- rowSums(kaiju_counts[,10:(9+sample_num)])
}else{
kaiju_counts$sum <- rowSums(kaiju_counts[,9:(8+sample_num)])
}

drops=c("sum")
kaiju_counts1=kaiju_counts[kaiju_counts$sum>=1,!(names(kaiju_counts) %in% drops)]
write.table(kaiju_counts1,file=output,sep = "\t",row.names =FALSE, col.names =FALSE)
