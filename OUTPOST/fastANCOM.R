# install.packages("kableExtra", dependencies = TRUE)
# library(kableExtra)
# https://rpubs.com/RChao/820127
library(fastANCOM)
library(ggplot2)
library(ggpubr)

# settings
args = commandArgs(T)
dp = args[1]
wkdir = args[2]
group1 = args[3]
group2 = args[4]
group1_index = args[5]
group2_index = args[6]
level = args[7]

data = as.matrix(read.csv(dp,row.names = 1,sep=",",check.names=FALSE))
group1_index = as.numeric(strsplit(group1_index, split = ",")[[1]])+1
group2_index = as.numeric(strsplit(group2_index, split = ",")[[1]])+1
data = data[c(group1_index,group2_index),]
group = c(rep(0,times = length(group1_index)),rep(1,times = length(group2_index)))

setwd(wkdir)

# fit model
fit <- fastANCOM(Y=data, x=group)
final_fit <- fit$results$final
write.table(final_fit,file = paste("ancom_biomarker.",group1,"_vs_",group2,".at_",level,".tsv",sep=""),sep="\t",col.names=NA)

#Results with plot
pdf(file=paste("ancom_biomarkers.dotplot.",group1,"_vs_",group2,".at_",level,".pdf",sep=""))
wp <- ggplot(data=final_fit, aes(x=1:dim(data)[2],y=Reject.number, color=REJECT, shape=REJECT)) + 
  geom_bar(aes(y=Reject.number), stat = "identity", width = 0.005, color='lightgrey') +
  # geom_vline(xintercept = 1:100, linetype='dashed', size=0.1) +
  geom_point(size=2) + theme_linedraw() + 
  labs(x='Microbe relative abundance rank',y='Reject number', title='ANCOM identified biomarkers') +
  scale_colour_brewer(palette = "Set1", direction=-1) + 
  theme(legend.position = c(0.95, 0.95), 
        legend.background = element_rect(fill = "white", colour = "black"))
print(wp)
dev.off()

pdf(file=paste("ancom_biomarkers.volcano.",group1,"_vs_",group2,".at_",level,".pdf",sep=""))
vp <-ggplot(data=final_fit,aes(x=log2FC,y=-log10(log2FC.pval), color=REJECT, shape=REJECT)) +
  geom_hline(yintercept=-log10(0.05/nrow(final_fit)),
             linetype="dashed",color="black",size=0.5)+ 
  geom_text(aes(label=ifelse(final_fit$log2FC.pval < 0.05, row.names(final_fit), NA)), vjust=1, nudge_x=0.0, nudge_y=0.0, size=3, na.rm=TRUE) +
  geom_point(size=2) + theme_linedraw() + labs(y='log10(p-value)', title='Volcano plot') +
  theme(legend.position='none') +
  scale_colour_brewer(palette = "Set1",direction=-1) +
  theme(legend.position = c(0.95, 0.95), 
        legend.background = element_rect(fill = "white", colour = "black"))
print(vp)
dev.off()




