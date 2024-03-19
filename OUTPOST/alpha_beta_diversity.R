library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)
library(wesanderson)
library(ggpubr)

# settings
args = commandArgs(T)
dp = args[1]
wkdir = args[2]
index_dir = args[3]
group_dir = args[4]
groups = args[5]
level = args[6]

# dp = r"{D:\CurrentProjects\OUTPOST\horsedonkey\alpha_beta_analysis\alpha_beta_horsedonkey_vs_hinny\horsedonkey.taxa_counts.rel_abun.genus.rmU.horsedonkey_vs_hinny.csv}"
# wkdir=r"{D:\CurrentProjects\OUTPOST\horsedonkey\alpha_beta_analysis\alpha_beta_horsedonkey_vs_hinny}"
# index_dir=r"{D:\CurrentProjects\OUTPOST\horsedonkey\alpha_beta_analysis\alpha_beta_horsedonkey_vs_hinny\index.txt}"
# group_dir=r"{D:\CurrentProjects\OUTPOST\horsedonkey\alpha_beta_analysis\alpha_beta_horsedonkey_vs_hinny\group.csv}"
# groups = "horsedonkey,hinny"


my_comparisons = strsplit(groups,",")
index=read.csv(index_dir,header = FALSE)
group=read.csv(group_dir,row.names = 1,sep=",")
spe = read.csv(dp,sep=",",check.names=FALSE)
spe = spe[ , !(names(spe) %in% c("X",""))]
rownames(spe)=index$V1

# must read files first then set working dir
setwd(wkdir)
################################################################################
# alpha_diversity
################################################################################
Shannon = diversity(spe)
Inv_Simpson = diversity(spe,"invsimpson")
Simpson = diversity(spe,"simpson")
taxa_num = specnumber(spe)
Pielou_evenness = Shannon/log(taxa_num)
Simpson_evenness = Simpson/taxa_num

# output
report = cbind(Shannon,Simpson, Inv_Simpson,Pielou_evenness, Simpson_evenness)
report = as.data.frame(report) %>% rownames_to_column(var = "ID")
report_dp = paste("alpha_diversity.at_",level,".csv",sep = "")
write.table(report,file = report_dp,sep=",",row.names = FALSE)

# draw alpha_diversity box plot
alpha_diversity = read.csv(report_dp,sep = ",")  ### change the input_table
colnames(alpha_diversity)[1] = "Group"
alpha_diversity$Group = group$Group
alpha_diversity = Filter(function(x)!all(is.na(x)), alpha_diversity)

alphas = colnames(alpha_diversity)
alphas = alphas[alphas!="Group"]

for (alpha in alphas){
  pdf(file=paste(alpha,"_alpha_diversity.at_",level,".pdf",sep=""))
  #draw scatter-plots and mark outliers
  p <- ggboxplot(alpha_diversity,x="Group",y=alpha,add="dotplot",color = "Group", palette = "jco") +
      theme(legend.position="right") +
      labs(title=alpha, y = "alpha_diversity") +
      theme_classic() +
      stat_compare_means(method = "wilcox",paired = F,comparisons = my_comparisons)
  print(p)
  dev.off()
}
################################################################################
# beta_diversity
################################################################################
bray.dist = vegdist(spe,method="bray",na.rm = TRUE)
bray.mat = as.matrix(bray.dist)
bray.dist.ln = vegdist(log1p(spe),na.rm = TRUE)
bray.mat.ln=as.matrix(bray.dist.ln)
spe.jaccard<-vegdist(spe,method="jaccard")
spe.jaccard.mat=as.matrix(spe.jaccard)

write.table(bray.mat,file= paste("bray.at_",level,".csv",sep = ""),sep=",",col.names=NA)
write.table(bray.mat.ln,file=paste("bray_ln.at_",level,".csv",sep = ""),sep=",",col.names=NA)
write.table(spe.jaccard.mat,file=paste("jaccard.at_",level,".csv",sep = ""),sep=",",col.names=NA)

#plot PCoA
beta_pcoa=function (dis_mat, metadata, groupID = "Group", ellipse = T,
          label = F, PCo = 12)
{
  p_list = c("ggplot2", "vegan", "ggrepel")
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    suppressWarnings(suppressMessages(library(p, character.only = T)))
  }
  idx = rownames(metadata) %in% rownames(dis_mat)
  metadata = metadata[idx, , drop = F]
  dis_mat = dis_mat[rownames(metadata), rownames(metadata)]
  sampFile = as.data.frame(metadata[, groupID], row.names = row.names(metadata))
  pcoa = cmdscale(dis_mat, k = 3, eig = T)
  points = as.data.frame(pcoa$points)
  eig = pcoa$eig
  sum_eig=sum(abs(eig))
  points = cbind(points, sampFile[rownames(points), ])
  colnames(points) = c("x", "y", "z", "group")
  if (PCo == 12) {
    p = ggplot(points, aes(x = x, y = y, color = group)) +
      labs(x = paste("PCo 1 (", format(100 * eig[1]/sum_eig,
      digits = 4), "%)", sep = ""), y = paste("PCo 2 (",
      format(100 * eig[2]/sum_eig, digits = 4), "%)",
      sep = ""), color = groupID)
  }
  if (PCo == 13) {
    p = ggplot(points, aes(x = x, y = z, color = group)) +
      labs(x = paste("PCo 1 (", format(100 * eig[1]/sum_eig,
       digits = 4), "%)", sep = ""), y = paste("PCo 3 (",
        format(100 * eig[3]/sum_eig, digits = 4), "%)",
       sep = ""), color = groupID)
  }
  if (PCo == 23) {
    p = ggplot(points, aes(x = y, y = z, color = group)) +
      labs(x = paste("PCo 2 (", format(100 * eig[2]/sum_eig,
       digits = 4), "%)", sep = ""), y = paste("PCo 3 (",
       format(100 * eig[3]/sum_eig, digits = 4), "%)",
      sep = ""), color = groupID)
  }

  p = p + geom_point(alpha = 0.7, size = 2) + theme_classic() +
    theme(text = element_text(family = "sans", size = 7))

  if (ellipse == T) {
    p = p + stat_ellipse(level = 0.95)
  }

  if (label == T) {
    p = p + geom_text_repel(label = paste(rownames(points)),
                            colour = "black", size = 1)
  }
  if (PCo == 123) {
    # library("scatterplot3d") # load
    # colors=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99")
    # colors <- colors[as.numeric(as.factor(points$group))]
    # p = scatterplot3d(points[,1:3],
    #                   main="3D-PCoA plot",
    #                   xlab = "PCo 1",
    #                   ylab = "PCo 2",
    #                   zlab = "PCo 3",
    #                   color=colors, pch = 16)
    library("car")
    # 3D plot with the regression plane
    p=scatter3d(x=points[,1],y=points[,2],z=points[,3],
              groups = as.factor(points$group),
              ellipsoid = TRUE)
    }
  return(list(one=p,two=pcoa))
}
# jaccard
p12=beta_pcoa(dis_mat=spe.jaccard.mat,metadata=group, groupID="Group", ellipse = T,label = T, PCo = 12)
p13=beta_pcoa(dis_mat=spe.jaccard.mat,metadata=group, groupID="Group", ellipse = T,label = T, PCo = 13)
p23=beta_pcoa(dis_mat=spe.jaccard.mat,metadata=group, groupID="Group", ellipse = T,label = T, PCo = 23)
ggsave(paste0(paste("PCoA12.jaccard.at_",level,".pdf",sep = "")), p12$one, width=89, height=56, units="mm")
ggsave(paste0(paste("PCoA13.jaccard.at_",level,".pdf",sep = "")), p13$one, width=89, height=56, units="mm")
ggsave(paste0(paste("PCoA23.jaccard.at_",level,".pdf",sep = "")), p23$one, width=89, height=56, units="mm")
# bray curtis
p12=beta_pcoa(dis_mat=bray.mat,metadata=group, groupID="Group", ellipse = T, label = T, PCo = 12)
p13=beta_pcoa(dis_mat=bray.mat,metadata=group, groupID="Group", ellipse = T, label = T, PCo = 13)
p23=beta_pcoa(dis_mat=bray.mat,metadata=group, groupID="Group", ellipse = T, label = T, PCo = 23)
ggsave(paste0(paste("PCoA12.bray.at_",level,".pdf",sep = "")), p12$one, width=89, height=56, units="mm")
ggsave(paste0(paste("PCoA13.bray.at_",level,".pdf",sep = "")), p13$one, width=89, height=56, units="mm")
ggsave(paste0(paste("PCoA23.bray.at_",level,".pdf",sep = "")), p23$one, width=89, height=56, units="mm")

# add P-value
beta_pcoa_stat=function (dis_mat, metadata, groupID = "Group", result = paste("beta_pcoa_stat_at",level,".txt",sep = ""),method="bray")
{
  p_list = c("vegan")
  for (p in p_list) {
    if (!requireNamespace(p)) {
      install.packages(p)
    }
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
  }
  metadata$group = metadata[[groupID]]
  # write.table(date(), file = result, append = T, sep = "\t",
  #             quote = F, row.names = F, col.names = F)
  idx = rownames(metadata) %in% rownames(dis_mat)
  metadata = metadata[idx, , drop = F]
  dis_mat = dis_mat[rownames(metadata), rownames(metadata)]

  da_adonis = function(sampleV) {
    sampleA = as.matrix(sampleV$sampA)
    sampleB = as.matrix(sampleV$sampB)
    design2 = subset(metadata, group %in% c(sampleA, sampleB))
    if (length(unique(design2$group)) > 1) {
      sub_dis_table = dis_table[rownames(design2), rownames(design2)]
      # the following sentence is the core
      # remove the diagnoise and remain the tri-angle
      sub_dis_table = as.dist(sub_dis_table, diag = FALSE,upper = FALSE)
      adonis_table = adonis(sub_dis_table ~ group, data = design2,permutations = 10000)
      adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
      adonis_pvalue = paste(sampleA, sampleB, adonis_pvalue,sep = ",")
      write.table(adonis_pvalue, file = result, append = TRUE, sep = ",", quote = F, row.names = F, col.names = F)
    }
  }

  dis_table = as.matrix(dis_mat)

  compare_data = as.vector(unique(metadata$group))
  len_compare_data = length(compare_data)
  for (i in 1:(len_compare_data - 1)) {
    for (j in (i + 1):len_compare_data) {
      tmp_compare = as.data.frame(cbind(sampA = compare_data[i],sampB = compare_data[j]))
      da_adonis(tmp_compare)
    }
  }

}

# use adonis to detect pair-wise difference and save
beta_pcoa_stat(spe.jaccard.mat, group, "Group", paste("beta_pcoa_P-value.jaccard.at_",level,".csv",sep = ""),"jaccard")
beta_pcoa_stat(bray.mat, group, "Group",paste("beta_pcoa_P-value.bray.at_",level,".csv",sep = ""),"bray")


# 3d pca
# library(pca3d)
# data( metabo )
# pca <- prcomp( metabo[,-1], scale.= TRUE )
# pca2d( pca, group= metabo[,1] )
# pca3d( pca,show.ellipses = T,group= metabo[,1] )
#
# pcaa=prcomp(spe.jaccard.mat,scale. = T)
# pca2d(spe.jaccard.mat,group = group$Group)
# pca3d( pcaa,show.ellipses = T,group= group$Group )
#
# p12=beta_pcoa(dis_mat=spe.jaccard.mat,metadata=group, groupID="Group", ellipse = T,
#               label = F, PCo = 12)
