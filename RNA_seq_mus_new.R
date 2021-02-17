##DEseq2差异基因分组分析
#分别整理各组分析表
expression_all<-read.delim('Expression_all.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                           check.names = F)
##分析D10NTC与饲喂d10，d20，d30，d35差异基因情况
#摘取不同组counts数据
D10NTC_D10TOA_count<-expression_all[,c(2,4,6,14,16,18)]
D20NTC_D20TOA_count<-expression_all[,c(20,22,24,32,34,36)]
D30NTC_D30TOA_count<-expression_all[,c(38,40,42,50,52,54)]
D30NTC_D35TOA_count<-expression_all[,c(38,40,42,56,58,60)]
###建立Deseq2流程的function函数
Deseq_method<-function(data,n,group1,group2){
  a<-factor(c(rep("group1",n),rep("group2",n)), levels = c("group1","group2"))
  b<-data.frame(row.names = colnames(data),a)
  require(DESeq2)
  dds<- DESeqDataSetFromMatrix(data, b, design= ~ a)
  dds<-dds[rowSums(counts(dds))>1,]
  dds<-DESeq(dds)
  res<- na.omit(results(dds, contrast=c("a", "group2", "group1")))
  res_1<-as.data.frame(res)
  return(res_1)
}
D10NTC_D10TOA_DEG<-as.data.frame(Deseq_method(D10NTC_D10TOA_count,3,D10NTC,D10TOA))
D20NTC_D20TOA_DEG<-as.data.frame(Deseq_method(D20NTC_D20TOA_count,3,D20NTC,D20TOA))
D30NTC_D30TOA_DEG<-as.data.frame(Deseq_method(D30NTC_D30TOA_count,3,D30NTC,D30TOA))
D30NTC_D35TOA_DEG<-as.data.frame(Deseq_method(D30NTC_D35TOA_count,3,D30NTC,D35TOA))
#添加上升下降趋势
D10NTC_D10TOA_DEG['significant']<-ifelse(
  D10NTC_D10TOA_DEG$pvalue < 0.05 & D10NTC_D10TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10NTC_D10TOA_DEG$pvalue < 0.05 & D10NTC_D10TOA_DEG$log2FoldChange < -1,"down","no"))
D20NTC_D20TOA_DEG['significant']<-ifelse(
  D20NTC_D20TOA_DEG$pvalue < 0.05 & D20NTC_D20TOA_DEG$log2FoldChange >1,"up",
  ifelse(D20NTC_D20TOA_DEG$pvalue < 0.05 & D20NTC_D20TOA_DEG$log2FoldChange < -1,"down","no"))
D30NTC_D30TOA_DEG['significant']<-ifelse(
  D30NTC_D30TOA_DEG$pvalue < 0.05 & D30NTC_D30TOA_DEG$log2FoldChange >1,"up",
  ifelse(D30NTC_D30TOA_DEG$pvalue < 0.05 & D30NTC_D30TOA_DEG$log2FoldChange < -1,"down","no"))
D30NTC_D35TOA_DEG['significant']<-ifelse(
  D30NTC_D35TOA_DEG$pvalue < 0.05 & D30NTC_D35TOA_DEG$log2FoldChange >1,"up",
  ifelse(D30NTC_D35TOA_DEG$pvalue < 0.05 & D30NTC_D35TOA_DEG$log2FoldChange < -1,"down","no"))
#ggplot2火山图可视化
#火山图绘图函数
volcano_plot<-function(data){
  data['significant']<-ifelse(data$pvalue < 0.05 & data$log2FoldChange >1,"up",
                              ifelse(data$pvalue < 0.05 & data$log2FoldChange < -1,"down","no"))
  require(ggplot2)
  plot<-ggplot(data,aes(x=log2FoldChange,y=-1*log10(pvalue)))+
    geom_point(aes(color=significant),size=2)+
    xlim(-5,5)+ylim(0,7)+labs(title="Volcano Plot",x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
    scale_color_manual(values =c("blue",'gray','red'))+
    geom_hline(yintercept=1.3,linetype=2,col="black")+
    geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
    theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
          axis.text.x = element_text(size=15,colour = 'black'),
          legend.text = element_text(size=15,colour = 'black'),legend.title = element_text(size=15,colour = 'black'))
  return(plot)
}
D10NTCvsD10TOA_plot<-volcano_plot(D10NTC_D10TOA_DEG)
D20NTCvsD20TOA_plot<-volcano_plot(D20NTC_D20TOA_DEG)
D30NTCvsD30TOA_plot<-volcano_plot(D30NTC_D30TOA_DEG)
D30NTCvsD35TOA_plot<-volcano_plot(D30NTC_D35TOA_DEG)
D10NTCvsD10TOA_plot
library(gridExtra)
all_plot<-grid.arrange(D10NTCvsD10TOA_plot,D20NTCvsD20TOA_plot,
             D30NTCvsD30TOA_plot,D30NTCvsD35TOA_plot,nrow=2,ncol=2)
ggsave("all_volcano plot.png",all_plot,width=10,height=12)
#分别摘取上调与下调基因数据框
D10NTCvsD10TOA_up<-D10NTC_D10TOA_DEG[which(D10NTC_D10TOA_DEG$significant == 'up'),]
D20NTCvsD20TOA_up<-D20NTC_D20TOA_DEG[which(D20NTC_D20TOA_DEG$significant == 'up'),]
D30NTCvsD30TOA_up<-D30NTC_D30TOA_DEG[which(D30NTC_D30TOA_DEG$significant == 'up'),]
D30NTCvsD35TOA_up<-D30NTC_D35TOA_DEG[which(D30NTC_D35TOA_DEG$significant == 'up'),]
D10NTCvsD10TOA_down<-D10NTC_D10TOA_DEG[which(D10NTC_D10TOA_DEG$significant == 'down'),]
D20NTCvsD20TOA_down<-D20NTC_D20TOA_DEG[which(D20NTC_D20TOA_DEG$significant == 'down'),]
D30NTCvsD30TOA_down<-D30NTC_D30TOA_DEG[which(D30NTC_D30TOA_DEG$significant == 'down'),]
D30NTCvsD35TOA_down<-D30NTC_D35TOA_DEG[which(D30NTC_D35TOA_DEG$significant == 'down'),]
#柱状图统计
detail<-read.delim('detail.txt')
library(reshape2)
detail_1<-melt(detail,variable.name = 'Type',value.name = 'Count')
library(ggplot2)
diff.stat.plot<-ggplot(detail_1,aes(x=Group,y=Count,fill=Type))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+coord_flip()+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=14,colour = 'black'),
        legend.title =element_text(size=14,colour = 'black'))
diff.stat.plot
ggsave('diff_stat.png',diff.stat.plot,width = 8,height = 8)
#pick the all DEGs
expression_all['Gene_ID']<-row.names(expression_all)
gene_detail<-expression_all[,c(67,77)]
D10NTCvsD10TOA_up['Gene_ID']<-row.names(D10NTCvsD10TOA_up)
D10NTCvsD10TOA_down['Gene_ID']<-row.names(D10NTCvsD10TOA_down)
D20NTCvsD20TOA_up['Gene_ID']<-row.names(D20NTCvsD20TOA_up)
D20NTCvsD20TOA_down['Gene_ID']<-row.names(D20NTCvsD20TOA_down)
D30NTCvsD30TOA_up['Gene_ID']<-row.names(D30NTCvsD30TOA_up)
D30NTCvsD30TOA_down['Gene_ID']<-row.names(D30NTCvsD30TOA_down)
D30NTCvsD35TOA_up['Gene_ID']<-row.names(D30NTCvsD35TOA_up)
D30NTCvsD35TOA_down['Gene_ID']<-row.names(D30NTCvsD35TOA_down)
up_gene_list<-merge(merge(merge(D10NTCvsD10TOA_up,D20NTCvsD20TOA_up,by = 'Gene_ID',all = T),
                    D30NTCvsD30TOA_up,by = 'Gene_ID',all = T),D30NTCvsD35TOA_up,by = 'Gene_ID',
                    all = T)
up_gene_list<-as.data.frame(up_gene_list[,1])
down_gene_list<-merge(merge(merge(D10NTCvsD10TOA_down,D20NTCvsD20TOA_up,by = 'Gene_ID',all = T),
                          D30NTCvsD30TOA_down,by = 'Gene_ID',all = T),D30NTCvsD35TOA_down,by = 'Gene_ID',
                    all = T)
down_gene_list<-as.data.frame(down_gene_list[,1])
up_gene_list$Gene_ID<-up_gene_list$`up_gene_list[, 1]`
down_gene_list$Gene_ID<-down_gene_list$`down_gene_list[, 1]`
all_diff_genes<-merge(up_gene_list,down_gene_list,by = 'Gene_ID',all = T)
#DEG heatmaps
expression_test_fpkm<-expression_all[,c(3,5,7,15,17,19,21,23,25,33,35,37,39,41,43,51,53,55,57,59,61)]
expression_test_fpkm['Gene_ID']<-row.names(expression_test_fpkm)
expression_diff<-merge(expression_test_fpkm,all_diff_genes,by = 'Gene_ID',all.y   = T)
row.names(expression_diff)<-expression_diff$Gene_ID
expression_diff<-expression_diff[,-c(1,23:24)]
names(expression_diff)<-c('D10NTC1','D10NTC2','D10NTC3','D10TOA1','D10TOA2','D10TOA3',
                          'D20NTC1','D20NTC2','D20NTC3','D20TOA1','D20TOA2','D20TOA3',
                          'D30NTC1','D30NTC2','D30NTC3','D30TOA1','D30TOA2','D30TOA3',
                          'D35TOA1','D35TOA2','D35TOA3')
library(pheatmap)
group<-read.delim('group.txt',row.names = 1)
col<-colorRampPalette(c("blue", "white", "red"))(256)
df_1<-as.matrix(scale(expression_diff))
diff_heatmap<-pheatmap(df_1,scale = 'row',color=col,show_rownames = F,cellwidth = 25,
                       annotation_col = group,treeheight_row = 30,cluster_cols = F)
diff_heatmap

##GO富集
##分别GO富集
library(clusterProfiler)
library(org.Mm.eg.db)
GO_D10NTC_D10TOA_DEG_up<- enrichGO(gene = row.names(D10NTCvsD10TOA_up),OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D10NTC_D10TOA_DEG_down<- enrichGO(gene = row.names(D10NTCvsD10TOA_down),OrgDb= org.Mm.eg.db,
                                     keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D20NTC_D20TOA_DEG_up<- enrichGO(gene = row.names(D20NTCvsD20TOA_up),OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D20NTC_D20TOA_DEG_down<- enrichGO(gene = row.names(D20NTCvsD20TOA_down),OrgDb= org.Mm.eg.db,
                                     keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D30NTC_D30TOA_DEG_up<- enrichGO(gene = row.names(D30NTCvsD30TOA_up),OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D30NTC_D30TOA_DEG_down<- enrichGO(gene = row.names(D30NTCvsD30TOA_down),OrgDb= org.Mm.eg.db,
                                     keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D30NTC_D35TOA_DEG_up<- enrichGO(gene = row.names(D30NTCvsD35TOA_up),OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D30NTC_D35TOA_DEG_down<- enrichGO(gene = row.names(D30NTCvsD35TOA_down),OrgDb= org.Mm.eg.db,
                                     keyType= 'ENSEMBL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.2)
GO_D10NTC_D10TOA_DEG_up_detail<-as.data.frame(GO_D10NTC_D10TOA_DEG_up)
GO_D10NTC_D10TOA_DEG_down_detail<-as.data.frame(GO_D10NTC_D10TOA_DEG_down)
GO_D20NTC_D20TOA_DEG_up_detail<-as.data.frame(GO_D20NTC_D20TOA_DEG_up)
GO_D20NTC_D20TOA_DEG_down_detail<-as.data.frame(GO_D20NTC_D20TOA_DEG_down)
GO_D30NTC_D30TOA_DEG_up_detail<-as.data.frame(GO_D30NTC_D30TOA_DEG_up)
GO_D30NTC_D30TOA_DEG_down_detail<-as.data.frame(GO_D30NTC_D30TOA_DEG_down)
GO_D30NTC_D35TOA_DEG_up_detail<-as.data.frame(GO_D30NTC_D35TOA_DEG_up)
GO_D30NTC_D35TOA_DEG_down_detail<-as.data.frame(GO_D30NTC_D35TOA_DEG_down)
#计算rich-factor并排序
#建立函数
richfactor<-function(data,n){
  data['total']<-n
  data['Richfactor']<- data$Count /data$total
  data_decrease<-data[order(data$Richfactor,decreasing = T),]
  return(data_decrease)
}
GO_D10NTC_D10TOA_DEG_up_detail<-richfactor(GO_D10NTC_D10TOA_DEG_up_detail,56)
GO_D10NTC_D10TOA_DEG_down_detail<-richfactor(GO_D10NTC_D10TOA_DEG_down_detail,115)
GO_D20NTC_D20TOA_DEG_up_detail<-richfactor(GO_D20NTC_D20TOA_DEG_up_detail,23)
GO_D20NTC_D20TOA_DEG_down_detail<-richfactor(GO_D20NTC_D20TOA_DEG_down_detail,9)
GO_D30NTC_D30TOA_DEG_up_detail<-richfactor(GO_D30NTC_D30TOA_DEG_up_detail,174)
GO_D30NTC_D30TOA_DEG_down_detail<-richfactor(GO_D30NTC_D30TOA_DEG_down_detail,76)
GO_D30NTC_D35TOA_DEG_up_detail<-richfactor(GO_D30NTC_D35TOA_DEG_up_detail,140)
GO_D30NTC_D35TOA_DEG_down_detail<-richfactor(GO_D30NTC_D35TOA_DEG_down_detail,130)
#GO气泡图汇总
#建立气泡图函数
bubbleplot<-function(data){
  require(ggplot2)
  plot<-ggplot(data ,aes(-1*log10(Richfactor),Description))+
    geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
          axis.text.x = element_text(size=12,colour = 'black'))
  return(plot)
}
GO_D10NTC_D10TOA_DEG_up_top20plot<-bubbleplot(GO_D10NTC_D10TOA_DEG_up_detail[c(1:20),])
GO_D10NTC_D10TOA_DEG_down_top20plot<-bubbleplot(GO_D10NTC_D10TOA_DEG_down_detail[c(1:20),])
GO_D20NTC_D20TOA_DEG_up_top20plot<-bubbleplot(GO_D20NTC_D20TOA_DEG_up_detail[c(1:20),])
GO_D20NTC_D20TOA_DEG_down_top20plot<-bubbleplot(GO_D20NTC_D20TOA_DEG_down_detail[c(1:20),])
GO_D30NTC_D30TOA_DEG_up_top20plot<-bubbleplot(GO_D30NTC_D30TOA_DEG_up_detail[c(1,3:21),])
GO_D30NTC_D30TOA_DEG_down_top20plot<-bubbleplot(GO_D30NTC_D30TOA_DEG_down_detail[c(1:20),])
GO_D30NTC_D35TOA_DEG_up_top20plot<-bubbleplot(GO_D30NTC_D35TOA_DEG_up_detail[c(1:2,4:21),])
GO_D30NTC_D35TOA_DEG_down_top20plot<-bubbleplot(GO_D30NTC_D35TOA_DEG_down_detail[c(1:20),])
ggsave('GO_D10NTC_D10TOA_DEG_up_top20plot.tiff',GO_D10NTC_D10TOA_DEG_up_top20plot,height = 6,width = 8)
ggsave('GO_D20NTC_D20TOA_DEG_up_top20plot.tiff',GO_D20NTC_D20TOA_DEG_up_top20plot,height = 6,width = 8)
ggsave('GO_D30NTC_D30TOA_DEG_up_top20plot.tiff',GO_D30NTC_D30TOA_DEG_up_top20plot,height = 6,width = 8)
ggsave('GO_D30NTC_D35TOA_DEG_up_top20plot.tiff',GO_D30NTC_D35TOA_DEG_up_top20plot,height = 6,width = 8)
ggsave('GO_D10NTC_D10TOA_DEG_down_top20plot.tiff',GO_D10NTC_D10TOA_DEG_down_top20plot,height = 6,width = 8)
ggsave('GO_D20NTC_D20TOA_DEG_down_top20plot.tiff',GO_D20NTC_D20TOA_DEG_down_top20plot,height = 6,width = 8)
ggsave('GO_D30NTC_D30TOA_DEG_down_top20plot.tiff',GO_D30NTC_D30TOA_DEG_down_top20plot,height = 6,width = 8)
ggsave('GO_D30NTC_D35TOA_DEG_down_top20plot.tiff',GO_D30NTC_D35TOA_DEG_down_top20plot,height = 6,width = 8)
library(gridGraphics)
all_bubbleplot_up<-grid.arrange(GO_D10NTC_D10TOA_DEG_up_top20plot,GO_D20NTC_D20TOA_DEG_up_top20plot,
                                GO_D30NTC_D30TOA_DEG_up_top20plot,GO_D30NTC_D35TOA_DEG_up_top20plot,nrow=2,ncol=2)
all_bubbleplot_up
write.csv(file = 'GO_D10NTC_D10TOA_DEG_up.csv',GO_D10NTC_D10TOA_DEG_up_detail)
write.csv(file = 'GO_D30NTC_D30TOA_DEG_up.csv',GO_D30NTC_D30TOA_DEG_up_detail)
#抗菌肽时间序列分析
anti<-read.delim('anti.txt')
library(clusterProfiler)
library(org.Mm.eg.db)
anti<-bitr(unique(anti$Gene_ID), fromType = "ENSEMBL",toType = c( "SYMBOL"),OrgDb = org.Mm.eg.db)
anti['Gene_ID']<-anti$ENSEMBL
anti_1<-merge(anti,expression_test_fpkm,by = 'Gene_ID',all.x = T)
anti_1<-anti_1[,-c(2,10:12,16:18)]
anti_1['D10NTC_mean']<-apply(anti_1[,c(3:5)], 1, mean)
anti_1['D10TOA_mean']<-apply(anti_1[,c(6:8)], 1, mean)
anti_1['D20TOA_mean']<-apply(anti_1[,c(9:11)], 1, mean)
anti_1['D30TOA_mean']<-apply(anti_1[,c(12:14)], 1, mean)
anti_1['D35TOA_mean']<-apply(anti_1[,c(15:17)], 1, mean)
anti_1['D10NTC']<-0
anti_1['D10TOA']<-log2(anti_1$D10TOA_mean/anti_1$D10NTC_mean)
anti_1['D20TOA']<-log2(anti_1$D20TOA_mean/anti_1$D10NTC_mean)
anti_1['D30TOA']<-log2(anti_1$D30TOA_mean/anti_1$D10NTC_mean)
anti_1['D35TOA']<-log2(anti_1$D35TOA_mean/anti_1$D10NTC_mean)
anti_2<-anti_1[-2,-c(1,3:22)]
anti_2<-anti_2[order(anti_2$D10TOA,decreasing = T),]
library(reshape2)
anti_2.1<-melt(anti_2,variable.name = 'Group',value.name = 'rate')
anti_2.1['label']<-ifelse(anti_2.1$Group == 'D35TOA',anti_2$SYMBOL,'')
library(ggplot2)
anti_2.1_plot<-ggplot(anti_2.1,aes(x=Group,y=rate,colour=SYMBOL,group=SYMBOL))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide gene plot",x='Group',y=expression(log[2](Fold-Change)))+geom_text(aes(label=label),size=3,vjust=1)+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text(size=13,colour = 'black'),
        legend.text =element_text(size=13,colour = 'black'))
anti_2.1_plot
ggsave('anti_2.1plot.png',anti_2.1_plot,width = 8,height = 8)
anti_3<-anti_2
row.names(anti_3)<-anti_3$SYMBOL
anti_3<-anti_3[,-1]
anti_3.1<-as.matrix(anti_3)
library(pheatmap)
anti_heatmap<-pheatmap(anti_3,cutree_rows = 6,cluster_cols = F,
                       fontsize_row = 15,fontsize_col = 15,angle_col = 0)
anti_heatmap
