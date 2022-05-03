library(grid) 
library(limma)
library(edgeR)
library(DESeq2)
library(pheatmap)

# settings
args = commandArgs(T)
dp = args[1]
wkdir = args[2]
index_dir = args[3]
data_type = args[4] # data_type = taxa or humann or function or species or pvalue
output = args[5]

# dp=r"{C:\CurrentProjects\GEMINI\horsedonkey.taxa_counts.rel_abun.taxaID.rmU.top30.fillmin.log10.csv}"
# wkdir=r"{C:\CurrentProjects\GEMINI}"
# index_dir=r"{C:\CurrentProjects\GEMINI\horsedonkey.taxa_counts.rel_abun.taxaID.rmU.top30.fillmin.log10.index}"
# data_type='taxa'
# output = 'heatmap.pdf'

data = read.csv(dp,sep = ',')
data = as.matrix(data[, 2:ncol(data)])
row_index=read.csv(index_dir,sep=',',header = FALSE)
rownames(data)=c(as.vector(as.matrix(row_index)))

setwd(wkdir)

draw_colnames_90 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_90 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))

## scale = column/row/none
# x=pheatmap(data,scale='none',fontsize_row = 6,cluster_rows=T, border_color = "NA",cluster_cols=T ,breaks=c(seq(0,0.05,by=((0.05-(0))/100))))
# x=pheatmap(data,scale='none',fontsize_row = 6,cluster_rows=F, border_color = "NA",cluster_cols=F,breaks=c(seq(-6,0,by=((0-(-6))/100))))
x=pheatmap(data,scale="none",fontsize_row = 6,cluster_rows=T, border_color = "NA",cluster_cols=T)
if(data_type=='species'){
  width=dim(data)[2]/6
  height= dim(data)[1]
}else if(data_type=='taxa'){
  width=dim(data)[2]/4
  height= dim(data)[1]/2
}else if(data_type=='humann'){
  width=dim(data)[2]/3
  height= dim(data)[1]*2
}else if(data_type=='function'){
  width=dim(data)[2]/8
  height= dim(data)[1]*0.7
}else if(data_type=='pvalue'){
  width=dim(data)[2]
  height= dim(data)[1]
}

save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(x, output ,width,height)
dev.off()

