library(reshape2)
library(ggplot2)

# settings
args = commandArgs(T)
dp = args[1]
output = args[2]

abun_df = read.table(dp,header = TRUE,sep=",",row.names = 1)

sample_v = as.character(rownames(abun_df))
abun_df$sample=sample_v

length=dim(abun_df)[2]-1
if (length==12){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928")
}else if(length>=16){
  aa=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B")
  cols=tolower(c(aa,aa,aa))
}else if(length==16){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23","#BDB76B")
}else if(length==15){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080","#6B8E23")
}else if(length==14){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000","#808080")
}else if(length==13){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493","#b15928","#000000")
}else if(length==11){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#FF1493")
}else if(length==10){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
}else if(length==9){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")
}else if(length==8){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00")
}else if(length==7){
  cols=c("#327CBD","#F52E1A","#FF8500","#A05CA7","#50A150","#B05B28","#E2E2E2")
}else if(length==6){
  cols=c("#327CBD","#F52E1A","#FF8500","#A05CA7","#50A150","#B05B28")
}else if(length==5){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")
}else if(length==4){
  cols=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")
}else if(length==3){
  cols=c("#a6cee3","#1f78b4","#b2df8a")
}else if(length==2){
  cols=c("#a6cee3","#1f78b4")
}else if(length==1){
  cols=c("#a6cee3")
}

df_long <- melt(abun_df, id.vars = "sample", variable.name = "taxa")

pdf(output,width = (length(sample_v) * 1),height = (length(sample_v)/2))
ggplot(df_long, aes(x = sample, y = value, fill = taxa))+
       geom_bar(stat = "identity",position = position_stack(reverse = TRUE))+
       theme(axis.text.x = element_text(angle = 45))+
       scale_fill_manual(values=cols)

dev.off()
