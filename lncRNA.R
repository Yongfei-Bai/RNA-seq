#lncRNA Heatmap
lncRNA<-read.delim('lncRNA stat.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                check.names = F)
group<-read.delim('design.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                  check.names = F)
lncRNA_1<-lncRNA[which(rowSums(lncRNA) > 0),]#去除全部为0的行
Blank<-apply(lncRNA_1[,1:3], 1,mean)
Bo18h<-apply(lncRNA_1[,4:6], 1,mean)
Bo24h<-apply(lncRNA_1[,7:9], 1,mean)
Bo33h<-apply(lncRNA_1[,10:12], 1,mean)
lncRNA_1.1<-data.frame(lncRNA_1,Blank,Bo18h,Bo24h,Bo33h)#求平均值之后形成矩阵
lncRNA_1.1<-lncRNA_1.1[,-c(1:12)]
df_1<-as.matrix(scale(lncRNA_1))
df_1.1<-as.matrix(scale(lncRNA_1.1))
library(gplots)
rc <- rainbow(nrow(df_1), start = 0, end = .3)
cc <- rainbow(ncol(df_1), start = 0, end = .3)
lncRNA_heatmap<-heatmap.2(df_1, scale = "row", col=bluered(100),trace = "none", 
          density.info=c("none"),RowSideColors=rc, ColSideColors=cc,
          keysize = 1,key.title = F,key.xlab = NULL,labRow = F)
lncRNA_heatmap
#
library(pheatmap)
lncRNA_heatmap_2<-pheatmap(df_1,scale = 'row',color=col,show_rownames = F,
                           cellwidth = 25,annotation_col = group,treeheight_row = 30)
lncRNA_heatmap_2.1<-pheatmap(df_1.1,scale = 'row',color=col,show_rownames = F,
                             cellwidth = 70,treeheight_row = 30)
lncRNA_heatmap_2
lncRNA_heatmap_2.1

#mRNA express diff heatmap
mRNA<-read.delim('mRNA_stat.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                   check.names = F)
mRNA_1<-mRNA[which(rowSums(mRNA) > 0),]#去除全部为0的行
Blank<-apply(mRNA_1[,1:3], 1,mean)
Bo18h<-apply(mRNA_1[,4:6], 1,mean)
Bo24h<-apply(mRNA_1[,7:9], 1,mean)
Bo33h<-apply(mRNA_1[,10:12], 1,mean)
mRNA_1.1<-data.frame(mRNA_1,Blank,Bo18h,Bo24h,Bo33h)#求平均值之后形成矩阵
mRNA_1.1<-mRNA_1.1[,-c(1:12)]
df_2<-as.matrix(scale(mRNA_1))
df_2.1<-as.matrix(scale(mRNA_1.1))
library(pheatmap)
mRNA_heatmap_2<-pheatmap(df_2,scale = 'row',color=col,show_rownames = F,
                         cellwidth = 25,annotation_col = group,treeheight_row = 30)
mRNA_heatmap_2.1<-pheatmap(df_2.1,scale = 'row',color=col,show_rownames = F,
                           cellwidth = 70,treeheight_row = 30)
mRNA_heatmap_2
mRNA_heatmap_2.1

#miRNA express diff heatmap
miRNA<-read.delim('mirna.exp_profile.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                 check.names = F)
miRNA_1<-miRNA[,-c(1:15)]
miRNA_1[miRNA_1 < 1]=NA #设置小于1的值为NA
miRNA_1<-na.omit(miRNA_1)#去除小于1的TPM值
Blank<-apply(miRNA_1[,1:3], 1,mean)
Bo18h<-apply(miRNA_1[,4:6], 1,mean)
Bo24h<-apply(miRNA_1[,7:9], 1,mean)
Bo33h<-apply(miRNA_1[,10:12], 1,mean)
miRNA_1.1<-data.frame(miRNA_1,Blank,Bo18h,Bo24h,Bo33h)
miRNA_1.1<-miRNA_1.1[,-c(1:12)]
df_3<-as.matrix(scale(miRNA_1))
df_3.1<-as.matrix(scale(miRNA_1.1))
library(pheatmap)
miRNA_heatmap_2<-pheatmap(df_3,scale = 'row',color=col,show_rownames = F,
                          cellwidth = 25,annotation_col = group,treeheight_row = 30)
miRNA_heatmap_2.1<-pheatmap(df_3.1,scale = 'row',color=col,show_rownames = F,
                            cellwidth = 70,treeheight_row = 30)
miRNA_heatmap_2
miRNA_heatmap_2.1
#GO与KO分析
library(BiocManager)
BiocManager::install("AnnotationHub")#已加载
BiocManager::install("org.Bt.eg.db")#加载牛属注释文件
BiocManager::install("clusterProfiler")
BiocManager::install("DO.db")
BiocManager::install("gridExtra")
BiocManager::install("igraph")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("GO.db")
BiocManager::install("GOSemSim")
BiocManager::install("qvalue")
BiocManager::install("rvcheck")
BiocManager::install("tidyr")
library(AnnotationHub)	#library导入需要使用的数据包
library(org.Bt.eg.db)   #牛属注释数据库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(GOstats)
gene<-read.delim('gene.txt', sep = '\t',stringsAsFactors = F, check.names = F)
gene<-read.csv('gene_1.csv')
goAnn <- get("org.Bt.egGO")
universe <- Lkeys(goAnn)
gene_1<-as.character(gene[,2])
ids <- bitr(gene_1,fromType = "SYMBOL",
            toType =c("ENTREZID",'ENSEMBL'),
            OrgDb = "org.Bt.eg.db")#转换accession num到ENTREZID
head(ids)
genes<- ids[,2]
head(genes)
ego2 <- enrichGO(gene         = ids$SYMBOL,
                 OrgDb         = org.Bt.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)
ego2<-as.data.frame(ego2)

ALL<- new("GOHyperGParams", geneIds=genes,universeGeneIds=universe,
          annotation="org.Bt.eg.db",
          ontology="BP",pvalueCutoff=0.05,conditional=FALSE,testDirection="over")
over <- hyperGTest(ALL)
library(Category)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {.sym <- merge(.ids, envir=org.Bt.egSYMBOL, 
                                                    ifnotfound=NA)
  + 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
  + 	paste(.sym, collapse=";")})
head(glist)
bp <- summary(over)



ekegg <- enrichKEGG(genes,keyType = "kegg",organism = "bta",
                    pvalueCutoff = 0.05,qvalueCutoff = 0.2,use_internal_data = FALSE)
ekegg_all<-as.data.frame(ekegg)

#
gene1<-read.delim('gene.txt', sep = '\t',stringsAsFactors = F, check.names = F)
KO<-read.delim('KO.txt', sep = '\t',stringsAsFactors = F, check.names = F)
mergetable<-merge(gene1,KO,by = 'gene',all.x = F)
mergetable1<-na.omit(mergetable)
write.csv(mergetable1,file = 'mergetable.csv')

#
gene_a<-read.delim('gene.txt', sep = '\t',stringsAsFactors = F, check.names = F)
gene_b<-read.delim('mRNA_1.txt', sep = '\t',stringsAsFactors = F, check.names = F)
gene_c<-read.delim('mRNA_2.txt', sep = '\t',stringsAsFactors = F, check.names = F)
gene_a1<-merge(gene_a,gene_b,by = 'GeneID',all.x = T)
gene_b1<-merge(gene_a1,gene_c,by = 'Symbol',all.x = T)
gene_b11<-gene_b1[-c(1:3770),]
geneb111<-merge(gene_b11,gene_a,by = 'GeneID',all = T)
geneb1111<-na.omit(geneb111)
write.csv(geneb1111,file = 'gene_1.csv')
gene_b<-gene_b[order(gene_b$Symbol,decreasing = F),]
gene_c<-gene_c[order(gene_c$Symbol,decreasing = F),]
gene_c['fa']<-c(1:21233)
gene_b_1<-gene_b[-c(1:11),]
gene_c_1<-gene_c[-c(1:755),-3]
mRNA_all<-merge(gene_b_1,gene_c_1,by = 'Symbol',all.x = T)


#
node1<-read.delim('Merged Network(1) default node.txt', sep = '\t',
                  stringsAsFactors = F, check.names = F)
node2<-read.delim('Merged Network(2) default node.txt',sep = '\t',
                  stringsAsFactors = F,check.names = F)
library(ggplot2)
BC_up_boxplot<-ggplot(node1,aes(x=type,y=BetweennessCentrality))+
  geom_boxplot(aes(fill=type))#BC_up指数
CC_up_boxplot<-ggplot(node1,aes(x=type,y=ClosenessCentrality))+
  geom_boxplot(aes(fill=type))#CC_up指数图
BC_down_boxplot<-ggplot(node2,aes(x=type,y=BetweennessCentrality))+
  geom_boxplot(aes(fill=type))#BC_down指数
CC_down_boxplot<-ggplot(node2,aes(x=type,y=ClosenessCentrality))+
  geom_boxplot(aes(fill=type))#CC_down指数图
BC_up_boxplot
CC_up_boxplot
BC_down_boxplot
CC_down_boxplot

#去除新组装lncRNA与预测miRNA重分析

lncRNA<-read.delim('lncRNA stat.txt',sep = '\t',stringsAsFactors = F, 
                   check.names = F)
lncRNA_diff<-read.delim('lncRNA.all.diff.txt',sep = '\t',stringsAsFactors = F, 
                        check.names = F)
group<-read.delim('design.txt', row.names = 1,sep = '\t',stringsAsFactors = F, 
                  check.names = F)
lncRNA_1<-merge(lncRNA,lncRNA_diff,by = 'id',all.x = T)
lncRNA_1<-na.omit(lncRNA_1)
row.names(lncRNA_1)<-lncRNA_1$id
lncRNA_1<-lncRNA_1[,-c(1,14)]
lncRNA_1<-lncRNA_1[which(rowSums(lncRNA_1) > 0),]#去除全部为0的行
df_1<-as.matrix(scale(lncRNA_1))
library(pheatmap)
library(gplots)
lncRNA_heatmap_2<-pheatmap(df_1,scale = 'row',color=col,show_rownames = F,
                           cellwidth = 25,annotation_col = group,treeheight_row = 30)
lncRNA_heatmap_2

#mRNA express diff heatmap
mRNA<-read.delim('mRNA_stat.txt',sep = '\t',stringsAsFactors = F, 
                 check.names = F)
mRNA_1<-mRNA[which(rowSums(mRNA) > 0),]#去除全部为0的行
mRNA_1<-mRNA_1[-c(16582:16774),]
df_2<-as.matrix(scale(mRNA_1))
library(pheatmap)
mRNA_heatmap_2<-pheatmap(df_2,scale = 'row',color=col,show_rownames = F,
                         cellwidth = 25,annotation_col = group,treeheight_row = 30)
mRNA_heatmap_2

#miRNA express diff heatmap
miRNA<-read.delim('mirna.exp_profile.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                  check.names = F)
miRNA_diff<-read.delim('miRNA.all.diff.txt',sep = '\t',stringsAsFactors = F, 
                       check.names = F)
miRNA_1<-miRNA[,-c(1:15)]
miRNA_1[miRNA_1 < 1]=NA #设置小于1的值为NA
miRNA_1<-na.omit(miRNA_1)#去除小于1的TPM值
miRNA_1['id']<-row.names(miRNA_1)
miRNA_2<-merge(miRNA_1,miRNA_diff,by = 'id',all.x = T)
miRNA_2<-na.omit(miRNA_2)
row.names(miRNA_2)<-miRNA_2$id
miRNA_2<-miRNA_2[,-c(1,14)]
head(miRNA_2)
df_3<-as.matrix(scale(miRNA_2))
library(pheatmap)
miRNA_heatmap_2<-pheatmap(df_3,scale = 'row',color=col,show_rownames = F,
                          cellwidth = 25,annotation_col = group,treeheight_row = 30)
miRNA_heatmap_2

#
node1<-read.delim('Merged Network(3) default node.txt', sep = '\t',
                  stringsAsFactors = F, check.names = F)
node2<-read.delim('Merged Network(4) default node.txt',sep = '\t',
                  stringsAsFactors = F,check.names = F)
library(ggplot2)
#BC_up指数
BC_up_boxplot<-ggplot(node1,aes(x=type,y=BetweennessCentrality))+
  geom_boxplot(aes(fill=type))+theme(legend.position = 'none',
                                    axis.title.y = element_text(size = 15, face = 'bold'),
                                    axis.text.x = element_text(size = 15, face = 'bold'))
#CC_up指数图
CC_up_boxplot<-ggplot(node1,aes(x=type,y=ClosenessCentrality))+
  geom_boxplot(aes(fill=type))+theme(legend.position = 'none',
                                     axis.title.y = element_text(size = 15, face = 'bold'),
                                     axis.text.x = element_text(size = 15, face = 'bold'))
#BC_down指数
BC_down_boxplot<-ggplot(node2,aes(x=type,y=BetweennessCentrality))+
  geom_boxplot(aes(fill=type))+theme(legend.position = 'none',
                                     axis.title.y = element_text(size = 15, face = 'bold'),
                                     axis.text.x = element_text(size = 15, face = 'bold'))
#CC_down指数图
CC_down_boxplot<-ggplot(node2,aes(x=type,y=ClosenessCentrality))+
  geom_boxplot(aes(fill=type))+theme(legend.position = 'none',
                                     axis.title.y = element_text(size = 15, face = 'bold'),
                                     axis.text.x = element_text(size = 15, face = 'bold'))
BC_up_boxplot
CC_up_boxplot
BC_down_boxplot
CC_down_boxplot
#
all.stat<-read.delim('all.stat.txt',sep = '\t',
                  stringsAsFactors = F,check.names = F)
miRNA.stat<-all.stat[,-c(4:7)]
mRNA.stat<-all.stat[,-c(2:3,6:7)]
lncRNA.stat<-all.stat[,-c(2:5)]
library(ggplot2)
library(reshape2)
miRNA.stat_1<-melt(miRNA.stat,variable.name ='type',value.name = 'Count')
miRNA.stat.plot<-ggplot(miRNA.stat_1,aes(x=group,y=Count,fill=type))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+coord_flip()
mRNA.stat_1<-melt(mRNA.stat,variable.name ='type',value.name = 'Count')
mRNA.stat.plot<-ggplot(mRNA.stat_1,aes(x=group,y=Count,fill=type))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+coord_flip()
lncRNA.stat_1<-melt(lncRNA.stat,variable.name ='type',value.name = 'Count')
lncRNA.stat.plot<-ggplot(lncRNA.stat_1,aes(x=group,y=Count,fill=type))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+coord_flip()
miRNA.stat.plot
mRNA.stat.plot
lncRNA.stat.plot