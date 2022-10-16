library(ggplot2)
library(wesanderson)
library(ggpubr)
library(Ipaper)

# settings
args = commandArgs(T)
dp = args[1]
output = args[2]
groups = args[3]
group_pair = args[4]

input_csv = read.csv(dp,sep=",")
colnames(input_csv)[1] = "treatment"
groups = sapply(groups, function(x) strsplit(x, ",")[[1]], USE.NAMES=FALSE)
input_csv$treatment = groups
input_csv = Filter(function(x)!all(is.na(x)), input_csv)

my_comparisons = strsplit(group_pair,",")
taxa = colnames(input_csv)
taxa = taxa[taxa!="treatment"]

for (taxo in taxa){
  pdf(file=paste(output,taxo,"boxplot.pdf",sep="."))
  #draw scatter-plots and mark outliers
  p<-#ggplot(input_csv, aes_string(x="treatment", y=taxo,fill="treatment")) +
    ggboxplot(input_csv,x="treatment",y=taxo,add="dotplot", color = "treatment", palette = "jco") +
    #geom_boxplot(outlier.colour=NULL, outlier.shape=NULL,outlier.size=NULL,coef=1.5)+
    theme(legend.position="right")+
    labs(title=taxo, y = "relative abundance")+
    #geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=1,colour="treatment")+
    theme_classic()+
    #the following compare is for increase/decrease comparison
    #stat_compare_means(method = "wilcox",paired = F, hide.ns = TRUE,method.args = list(alternative=alter))
    stat_compare_means(method = "wilcox",paired = F,comparisons = my_comparisons)
  print(p)
  dev.off()
}
try(dev.off())



