library(mixOmics)
library(gridExtra) # PCA sample plot with density
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)

#---------------------------------------------------------------------
# Principal component analysis (PCA) with density plots per component
#---------------------------------------------------------------------
Scatter_Density <- function(data = data, batch = batch, trt = NULL, expl.var = expl.var,
                            xlim = xlim, ylim = ylim, batch.legend.title = 'Batch', 
                            trt.legend.title = 'Treatment', density.lwd = 0.2,
                            title = NULL, title.cex = 1.5, legend.cex = 0.7, legend.title.cex =0.75){
  data = as.data.frame(data)
  batch = as.factor(batch)
  trt = as.factor(trt)
  if(nlevels(trt) >= 2){
    pMain <- ggplot(data = data, aes(x = data[ ,1], y = data[ ,2], colour = batch, shape = trt)) + 
      geom_point() + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
      ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
      scale_color_manual(values = color.mixo(1:10)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title, shape = trt.legend.title)
    
    pTop <- ggplot(data,aes(x = data[ ,1], fill = batch, linetype = trt)) + 
      geom_density(size = density.lwd, alpha = 0.5) + ylab('Density') + 
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
            plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
            axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(xlim[1], xlim[2]) + labs(title = title)
    
    pRight <- ggplot(data, aes(x=data[ ,2], fill = batch, linetype = trt)) + 
      geom_density(size = density.lwd,alpha = 0.5) +  coord_flip() + ylab('Density') +
      theme(axis.title.x = element_text(size = rel(0.8)), 
            axis.title.y = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(ylim[1], ylim[2])
    
  }else{
    pMain <- ggplot(data = data, aes(x = data[ ,1], y=data[ ,2], colour = batch)) + 
      geom_point(shape = 16) + xlab(paste0('PC1: ', round(as.numeric(expl.var[1])*100), '% expl.var')) + 
      ylab(paste0('PC2: ', round(as.numeric(expl.var[2])*100), '% expl.var')) + 
      scale_color_manual(values = color.mixo(1:10)) + theme_bw() + xlim(xlim[1], xlim[2]) + 
      ylim(ylim[1], ylim[2]) + labs(colour = batch.legend.title)
    
    pTop <- ggplot(data, aes(x = data[ ,1], fill = batch)) + 
      geom_density(size = density.lwd, alpha=0.5) + ylab('Density') + 
      theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.8)), 
            plot.title = element_text(hjust = 0.5, size = rel(title.cex)), 
            axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(xlim[1], xlim[2]) + labs(title = title)
    
    pRight <- ggplot(data, aes(x=data[ ,2], fill = batch)) + 
      geom_density(size = density.lwd, alpha = 0.5) +  coord_flip() + ylab('Density') +
      theme(axis.title.x = element_text(size = rel(0.8)), 
            axis.title.y = element_blank(), axis.line = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.background = element_blank()) + scale_fill_manual(values = color.mixo(1:10)) +
      xlim(ylim[1], ylim[2])
  }
  
  
  g <- ggplotGrob(pMain + theme(legend.position = 'right', legend.box = 'horizontal',
                                legend.direction = 'vertical', 
                                legend.key.height = unit(0.2, 'cm'),
                                legend.key.width = unit(0.1, 'cm'),
                                legend.title = element_text(size = rel(legend.title.cex)),
                                legend.spacing.x = unit(0.1, 'cm'),
                                legend.spacing.y = unit(0.1, 'cm'),
                                legend.text = element_text(size = rel(legend.cex))))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  
  grid.arrange(pTop + theme(legend.position = 'none'), legend, pMain + 
                 theme(legend.position = 'none'), pRight + theme(legend.position = 'none'), 
               ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  
}
percentile_norm = function(data = data, batch = batch, trt = trt){
  batch = as.factor(batch)
  trt = as.factor(trt)
  trt.list = list()
  data.pn.df = data.frame()
  for(i in 1:nlevels(batch)){
    trt.each.b = trt[batch == levels(batch)[i]]
    trt.list[[i]] = trt.each.b
    data.each.b.pn = percentileofscore(data[batch == levels(batch)[i],], 
                                       which(trt.each.b == levels(trt.each.b)[1]))
    data.pn.df = rbind(data.pn.df,data.each.b.pn)
  }
  names(trt.list) = levels(batch)
  data.pn.df.reorder = data.pn.df[rownames(data), ]
  return(data.pn.df.reorder)
}

# library(knitr)
# library(pheatmap) # heatmap
# library(xtable) # table
# library(vegan) # RDA
# library(AgiMicroRna) # RLE plot
# library(cluster) # silhouette coefficient
# library(variancePartition) # variance calculation
# library(pvca) # PVCA
# library(ruv) # RUVIII
# library(lmerTest) # lmer
# library(bapred) # FAbatch

# settings
args = commandArgs(T)
rm_batch_effect = args[1]
dp = args[2]
output = args[3]
config = args[4]
before_plot = args[5]
after_plot = args[6]
# rm_batch_effect = 'Combat'
# dp = 'asian_old_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.has_batch_effect.csv'
# dp = 'human64_3_batch_effect.taxa_counts.rel_abun.class.rmU.has_batch_effect.csv'
# output = 'asian_old_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.csv'
# config = 'OUTPOST_config_3_batch_effect.tsv'
# before_plot = 'asian_old_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.before_rm_batch_effect.pdf'
# after_plot = 'asian_old_genefamilies_uniref90names_relab_rxn_unstratified.named.rel_abun_format.after_rm_batch_effect.pdf'
# setwd(r'{C:\Users\Hui\Nutstore\.nutstore_eWloYW5nam9lQGZveG1haWwuY29t\CurrentProjects}')
# df = as.matrix(read.csv(r'{C:\Users\Hui\Nutstore\.nutstore_eWloYW5nam9lQGZveG1haWwuY29t\CurrentProjects\OUTPOST\batch_effect\allSamples_genefamilies_uniref90names_relab_ko_unstratified.named.rel_abun_format.csv}'
#                      ,row.names = 1))
# config_df = read.csv(r'{C:\Users\Hui\Nutstore\.nutstore_eWloYW5nam9lQGZveG1haWwuY29t\CurrentProjects\OUTPOST\test_OUTPOST\OUTPOST_config_3.tsv}'
#                   ,header=T,sep = '\t')
df = as.matrix(read.csv(dp, row.names = 1))
config_df = read.csv(config,header=T,sep = '\t')
batch = as.factor(config_df$batch)
names(batch) = as.vector(config_df$samples)
batch = batch[row.names(df)]

# 1.2.2添加偏移
# 缩放之后的数
# 我们需要向所有计数数据添 1 的偏移量，以处理 CLR 转换的零
# 根据尺度不变原理（Aitchison 1986），它通过 CLR 转换对原始计数或 TSS 数据返回相同的结
df = df + 0.000001

# 1.2.3中心对数比变
# 微生物组数据是组成的，并且具有不同的库大小。对此类数据使用标准统计方法可能会导致虚假结果，因此必须进一步转换数据。CLR 是选择的转换
df.clr = logratio.transfo(df, logratio = 'CLR')
class(df.clr) = 'matrix' 

# 2. detect batch-effect
# 2.1主成分分  (PCA) 与每个成分的密度
df.pca.before = pca(df.clr, ncomp = 3)

# heatmap
x_min = min(df.pca.before$variates$X[,'PC1']) * 1.2
x_max = max(df.pca.before$variates$X[,'PC1']) * 1.2
y_min = min(df.pca.before$variates$X[,'PC2']) * 1.2
y_max = max(df.pca.before$variates$X[,'PC2']) * 1.2

pdf(before_plot, width = 12, height = 12)
Scatter_Density(data = df.pca.before$variates$X, batch = batch, 
                expl.var = df.pca.before$prop_expl_var$X, 
                xlim = c(x_min,x_max), ylim = c(y_min,y_max), 
                batch.legend.title = 'batch', 
                title = 'Before batch effect correction')
dev.off()

if (rm_batch_effect == 'Combat'){
  df.after = try( t(ComBat(t(df.clr), batch = batch, par.prior = F, prior.plots = F)))
  if ('try-error' %in% class(df.after)){
    print("OUTPOST: rm_batch_effect: cannot use 'Combat', use 'Limma' instead, because")
    print(df.after)
    df.after = t(removeBatchEffect(t(df.clr), batch = batch))
    title_ = 'After Limma batch effect correction'
  }else{
  title_ = 'After Combat batch effect correction'
  }
}else{
  df.after = t(removeBatchEffect(t(df.clr), batch = batch))
  title_ = 'After Limma batch effect correction'
}

df.pca.after <- pca(df.after, ncomp = 3)

pdf(after_plot, width = 12, height = 12)
x_min = min(df.pca.after$variates$X[,'PC1']) * 1.2
x_max = max(df.pca.after$variates$X[,'PC1']) * 1.2
y_min = min(df.pca.after$variates$X[,'PC2']) * 1.2
y_max = max(df.pca.after$variates$X[,'PC2']) * 1.2
Scatter_Density(data = df.pca.after$variates$X, batch = batch,
                expl.var = df.pca.after$prop_expl_var$X,
                xlim = c(x_min,x_max), ylim = c(y_min,y_max),
                batch.legend.title = 'batch',
                title = title_)
dev.off()

# write tables
write.csv(df.after, output, row.names = T)





