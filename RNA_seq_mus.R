expression<-read.delim('Expression.txt')
diff<-read.delim('diff_gene_of_all_groups.txt')
expression<-expression[,-46]
diff<-diff[,-44]
expression_1<-merge(expression,diff,by = 'Gene_ID',all.x = T)
expression_1<-na.omit(expression_1)
write.csv(expression_1,file = 'Expression_1.csv')
expression_1.1<-expression_1[which(expression_1$D10NTC_vs_D10TOA_de == 'up'),]
anti<-read.delim('anti.txt')
expression_2<-expression[,-(32:45)]
expression_2['D10NTC_mean']<-apply(expression_2[,2:4], 1,mean)
expression_2['D10TOA_mean']<-apply(expression_2[,8:10], 1, mean)
expression_2['D20TOA_mean']<-apply(expression_2[,17:19], 1, mean)
expression_2['D30TOA_mean']<-apply(expression_2[,26:28], 1, mean)
expression_2['D35TOA_mean']<-apply(expression_2[,29:31], 1, mean)
expression_2<-expression_2[,-c(2:31)]
expression_2.1<-merge(expression_2,anti,by = 'Gene_ID',all.x = T)
expression_2.1<-na.omit(expression_2.1)
row.names(expression_2.1)<-expression_2.1$Name
write.csv(expression_2.1,file = 'expression_2.csv')
#各组存在差异的抗菌肽趋势分析
expression_2.1['D10NTC']<-0
expression_2.1['D10TOA']<-log2(expression_2.1$D10TOA_mean/expression_2.1$D10NTC_mean)
expression_2.1['D20TOA']<-log2(expression_2.1$D20TOA_mean/expression_2.1$D10NTC_mean)
expression_2.1['D30TOA']<-log2(expression_2.1$D30TOA_mean/expression_2.1$D10NTC_mean)
expression_2.1['D35TOA']<-log2(expression_2.1$D35TOA_mean/expression_2.1$D10NTC_mean)
expression_2.2<-expression_2.1[,c(37,39:43)]
library(reshape2)
expression_2.2.1<-melt(expression_2.2,variable.name = 'Group',value.name = 'rate')
library(ggplot2)
p<-ggplot(expression_2.2.1,aes(x=Group,y=rate,colour=Name,group=Name))+geom_line()
p
#正趋势4基因折线图
expression_2.2.2<-expression_2.1[c(6,8,9,11),c(37,39:43)]
library(reshape2)
expression_2.2.2<-melt(expression_2.2.2,variable.name = 'Group',value.name = 'rate')
library(ggplot2)
p2<-ggplot(expression_2.2.2,aes(x=Group,y=rate,colour=Name,group=Name))+geom_line()
p2

#
sus<-read.delim('select sus_1.txt')
mus<-read.delim('select mus_1.txt')
same_gene<-merge(sus,mus,by = intersect(names(sus)[1],names(mus)[1]),all.x = T)
same_gene<-same_gene[,-c(3:5)]
library(org.Mm.eg.db)
gene_1<-toTable(org.Mm.egACCNUM)
gene_2<-toTable(org.Mm.egSYMBOL)
expression_2.2['symbol']<-expression_2.2$Name
tran1<-merge(expression_2.2,gene_2, by='symbol',all.x=T)
tran1_1<-na.omit(tran1)
tran2<-merge(tran1_1,gene_1,by='gene_id',all.x=T)
write.csv(tran2,file = 'tran2.csv')

##DEseq2差异基因分组分析
#分别整理各组分析表
expression_all<-read.delim('Expression_all.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                           check.names = F)
##分析D10NTC与饲喂d10，d20，d30，d35差异基因情况
#摘取不同组counts数据
D10NTC_D20TOA_count<-expression_all[,c(2,4,6,32,34,36)]
D10NTC_D30TOA_count<-expression_all[,c(2,4,6,50,52,54)]
D10NTC_D35TOA_count<-expression_all[,c(2,4,6,56,58,60)]
#编写样本组lib
condition_1 <- factor(c(rep("D10NTC",3),rep("D20TOA",3)), levels = c("D10NTC","D20TOA"))#一定为factor格式
condition_2 <- factor(c(rep("D10NTC",3),rep("D30TOA",3)), levels = c("D10NTC","D30TOA"))
condition_3 <- factor(c(rep("D10NTC",3),rep("D35TOA",3)), levels = c("D10NTC","D35TOA"))
#结合形成样本组数据框
coldata_1<-data.frame(row.names = colnames(D10NTC_D20TOA_count),condition_1)
coldata_2<-data.frame(row.names = colnames(D10NTC_D30TOA_count),condition_2)
coldata_3<-data.frame(row.names = colnames(D10NTC_D35TOA_count),condition_3)
#dds=DESeqDataSet Object构建dss模型
library(DESeq2)
#D10NTC与D20TOA
dds_1<- DESeqDataSetFromMatrix(D10NTC_D20TOA_count, coldata_1, design= ~ condition_1)
dds_1 <- DESeq(dds_1)
#总体结果查看
res_1<- results(dds_1, contrast=c("condition_1", "D20TOA", "D10NTC"))
res_1<-na.omit(res_1)
write.csv(res_1,file="D10NTC_D20TOA_results.csv")
res_1_result<-as.data.frame(res_1)
#差异基因输出
D10NTC_D20TOA_DEG_up<-subset(res_1, pvalue < 0.05 & log2FoldChange >1)#DEG上调
D10NTC_D20TOA_DEG_down<-subset(res_1, pvalue < 0.05 & log2FoldChange < -1)#DEG下调
#D10NTC与D30TOA
dds_2<- DESeqDataSetFromMatrix(D10NTC_D30TOA_count, coldata_2, design= ~ condition_2)
dds_2 <- DESeq(dds_2)
#总体结果查看
res_2<- results(dds_2, contrast=c("condition_2", "D30TOA", "D10NTC"))
res_2<-na.omit(res_2)
#差异基因输出
D10NTC_D30TOA_DEG_up<-subset(res_2, pvalue < 0.05 & log2FoldChange > 1)#DEG上调
D10NTC_D30TOA_DEG_down<-subset(res_2, pvalue < 0.05 & log2FoldChange < -1)#DEG下调
#D10NTC与D35TOA
dds_3<- DESeqDataSetFromMatrix(D10NTC_D35TOA_count, coldata_3, design= ~ condition_3)
dds_3<- DESeq(dds_3)
#总体结果查看
res_3<-results(dds_3, contrast=c("condition_3", "D35TOA", "D10NTC"))
res_3<-na.omit(res_3)
#差异基因输出
D10NTC_D35TOA_DEG_up<-subset(res_3, pvalue < 0.05 & log2FoldChange >1)#DEG上调
D10NTC_D35TOA_DEG_down<-subset(res_3, pvalue < 0.05 & log2FoldChange< -1)#DEG下调
#火山图
library(readxl)
D10NTC_D10TOA_DEG<-read_xlsx('D10NTC_vs_D10TOA.DESeq.xlsx')
D10NTC_D20TOA_DEG<-as.data.frame(res_1)
D10NTC_D30TOA_DEG<-as.data.frame(res_2)
D10NTC_D35TOA_DEG<-as.data.frame(res_3)
#添加上升下降趋势
D10NTC_D10TOA_DEG<-D10NTC_D10TOA_DEG[,c(1,3:7)]
D10NTC_D10TOA_DEG<-D10NTC_D10TOA_DEG[-(which(D10NTC_D10TOA_DEG$baseMeanA == 0 | D10NTC_D10TOA_DEG$baseMeanB ==0)),]#删除带INF值
D10NTC_D10TOA_DEG$log2FoldChange<-as.numeric(D10NTC_D10TOA_DEG$log2FoldChange)
D10NTC_D10TOA_DEG['significant']<-ifelse(
  D10NTC_D10TOA_DEG$pval < 0.05 & D10NTC_D10TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10NTC_D10TOA_DEG$pval < 0.05 & D10NTC_D10TOA_DEG$log2FoldChange < -1,"down","no"))
D10NTC_D20TOA_DEG['significant']<-ifelse(
  D10NTC_D20TOA_DEG$pvalue < 0.05 & D10NTC_D20TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10NTC_D20TOA_DEG$pvalue < 0.05 & D10NTC_D20TOA_DEG$log2FoldChange < -1,"down","no"))
D10NTC_D30TOA_DEG['significant']<-ifelse(
  D10NTC_D30TOA_DEG$pvalue < 0.05 & D10NTC_D30TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10NTC_D30TOA_DEG$pvalue < 0.05 & D10NTC_D30TOA_DEG$log2FoldChange < -1,"down","no"))
D10NTC_D35TOA_DEG['significant']<-ifelse(
  D10NTC_D35TOA_DEG$pvalue < 0.05 & D10NTC_D35TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10NTC_D35TOA_DEG$pvalue < 0.05 & D10NTC_D35TOA_DEG$log2FoldChange < -1,"down","no"))
#ggplot2火山图可视化
library(ggplot2)
D10NTCvsD10TOA_plot<-ggplot(D10NTC_D10TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pval)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10NTC_vs_D10TOA Volcano Plot",
                              x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),
        legend.text = element_text(size=15,colour = 'black'),legend.title = element_text(size=15,colour = 'black'))
D10NTCvsD20TOA_plot<-ggplot(D10NTC_D20TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10NTC_vs_D20TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10NTCvsD30TOA_plot<-ggplot(D10NTC_D30TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10NTC_vs_D30TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10NTCvsD35TOA_plot<-ggplot(D10NTC_D35TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10NTC_vs_D35TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10NTCvsD10TOA_plot
D10NTCvsD20TOA_plot
D10NTCvsD30TOA_plot
D10NTCvsD35TOA_plot
ggsave("D10NTCvsD10TOA_volcano plot.png",D10NTCvsD10TOA_plot,width=8,height=8)
ggsave("D10NTCvsD20TOA_volcano plot.png",D10NTCvsD20TOA_plot,width=8,height=8)
ggsave("D10NTCvsD30TOA_volcano plot.png",D10NTCvsD30TOA_plot,width=8,height=8)
ggsave("D10NTCvsD35TOA_volcano plot.png",D10NTCvsD35TOA_plot,width=8,height=8)

#得到的上升下降的DEG转换为数据框
D10NTC_D20TOA_DEG_up<-as.data.frame(D10NTC_D20TOA_DEG_up)
D10NTC_D20TOA_DEG_down<-as.data.frame(D10NTC_D20TOA_DEG_down)
D10NTC_D30TOA_DEG_up<-as.data.frame(D10NTC_D30TOA_DEG_up)
D10NTC_D30TOA_DEG_down<-as.data.frame(D10NTC_D30TOA_DEG_down)
D10NTC_D35TOA_DEG_up<-as.data.frame(D10NTC_D35TOA_DEG_up)
D10NTC_D35TOA_DEG_down<-as.data.frame(D10NTC_D35TOA_DEG_down)
D10NTC_D20TOA_DEG_up['Gene_ID']<-row.names(D10NTC_D20TOA_DEG_up)
D10NTC_D20TOA_DEG_down['Gene_ID']<-row.names(D10NTC_D20TOA_DEG_down)
D10NTC_D30TOA_DEG_up['Gene_ID']<-row.names(D10NTC_D30TOA_DEG_up)
D10NTC_D30TOA_DEG_down['Gene_ID']<-row.names(D10NTC_D30TOA_DEG_down)
D10NTC_D35TOA_DEG_up['Gene_ID']<-row.names(D10NTC_D35TOA_DEG_up)
D10NTC_D35TOA_DEG_down['Gene_ID']<-row.names(D10NTC_D35TOA_DEG_down)
#获得symbol列表
expression_all['Gene_ID']<-row.names(expression_all)
gene_detail<-expression_all[,c(67:75,77)]
#merge函数结合数据框
D10NTC_D20TOA_DEG_up<-merge(D10NTC_D20TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10NTC_D20TOA_DEG_down<-merge(D10NTC_D20TOA_DEG_down,gene_detail,by='Gene_ID',all.x=T)
D10NTC_D30TOA_DEG_up<-merge(D10NTC_D30TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10NTC_D30TOA_DEG_down<-merge(D10NTC_D30TOA_DEG_down,gene_detail,by='Gene_ID',all.x=T)
D10NTC_D35TOA_DEG_up<-merge(D10NTC_D35TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10NTC_D35TOA_DEG_down<-merge(D10NTC_D35TOA_DEG_down,gene_detail,by='Gene_ID',all.x=T)
#输出数据csv文件
write.csv(D10NTC_D20TOA_DEG_up,file = 'D10NTC_D20TOA_DEG_up_gene.csv',
          row.names = D10NTC_D20TOA_DEG_up$Gene_ID)
write.csv(D10NTC_D20TOA_DEG_down,file = 'D10NTC_D20TOA_DEG_down_gene.csv',
          row.names = D10NTC_D20TOA_DEG_down$Gene_ID)
write.csv(D10NTC_D30TOA_DEG_up,file = 'D10NTC_D30TOA_DEG_up_gene.csv',
          row.names = D10NTC_D30TOA_DEG_up$Gene_ID)
write.csv(D10NTC_D30TOA_DEG_down,file = 'D10NTC_D30TOA_DEG_down_gene.csv',
          row.names = D10NTC_D30TOA_DEG_down$Gene_ID)
write.csv(D10NTC_D35TOA_DEG_up,file = 'D10NTC_D35TOA_DEG_up_gene.csv',
          row.names = D10NTC_D20TOA_DEG_up$Gene_ID)
write.csv(D10NTC_D35TOA_DEG_down,file = 'D10NTC_D35TOA_DEG_down_gene.csv',
          row.names = D10NTC_D35TOA_DEG_down$Gene_ID)
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
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=13,colour = 'black'))
diff.stat.plot
ggsave('diff_stat.png',diff.stat.plot,width = 8,height = 8)
#查找4组相同上升，下降基因(UpSetR)
library(UpSetR)
library(readxl)
D10NTC_D10TOA_DEG_up<-read_xlsx('D10NTC_vs_D10TOA.DESeq.Up.xlsx')
D10NTC_D10TOA_DEG_down<-read_xlsx('D10NTC_vs_D10TOA.DESeq.Down.xlsx')
listinput_up<-list(D10NTC_D10TOA_DEG_up= D10NTC_D10TOA_DEG_up$id,
                     D10NTC_D20TOA_DEG_up = D10NTC_D20TOA_DEG_up$Gene_ID,
                     D10NTC_D30TOA_DEG_up = D10NTC_D30TOA_DEG_up$Gene_ID,
                     D10NTC_D35TOA_DEG_up = D10NTC_D35TOA_DEG_up$Gene_ID)
listinput_down<-list(D10NTC_D10TOA_DEG_down =D10NTC_D10TOA_DEG_down$id,
                     D10NTC_D20TOA_DEG_down =D10NTC_D20TOA_DEG_down$Gene_ID,
                     D10NTC_D30TOA_DEG_down =D10NTC_D30TOA_DEG_down$Gene_ID,
                     D10NTC_D35TOA_DEG_down =D10NTC_D35TOA_DEG_down$Gene_ID)
listinput_up_plot <- upset(fromList(listinput_up),nsets = 5, order.by = "freq")
listinput_down_plot<-upset(fromList(listinput_down),nsets = 5, order.by = "freq")
listinput_up_plot
listinput_down_plot
##热图
#获得所有组上调基因EMBL编号
D10NTC_D10TOA_DEG_up['Gene_ID']<-D10NTC_D10TOA_DEG_up$id
D10NTC_D10TOA_DEG_down['Gene_ID']<-D10NTC_D10TOA_DEG_down$id
up_gene_list<-merge(D10NTC_D20TOA_DEG_up,D10NTC_D10TOA_DEG_up,by = 'Gene_ID',all = T )
up_gene_list<-merge(up_gene_list,D10NTC_D30TOA_DEG_up,by = 'Gene_ID',all = T)
up_gene_list<-merge(up_gene_list,D10NTC_D35TOA_DEG_up,by = 'Gene_ID',all = T)
row.names(up_gene_list)<-up_gene_list$Gene_ID
up_gene_list<-up_gene_list[,1]
up_gene_list<-as.data.frame(up_gene_list)
up_gene_list['Gene_ID']<-up_gene_list$up_gene_list
#获得所有组下调基因EMBL编号
down_gene_list<-merge(D10NTC_D20TOA_DEG_down,D10NTC_D10TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list<-merge(down_gene_list,D10NTC_D30TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list<-merge(down_gene_list,D10NTC_D35TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list<-down_gene_list[,1]
down_gene_list<-as.data.frame(down_gene_list)
down_gene_list['Gene_ID']<-down_gene_list$down_gene_list
all_gene_diff<-merge(up_gene_list,down_gene_list,by = 'Gene_ID',all = T)
#获得所有基因fpkm值
expression_test_fpkm<-expression_all[,c(3,5,7,15,17,19,33,35,37,51,53,55,57,59,61,77)]
#get diff gene fpkm value
expression_diff<-merge(expression_test_fpkm,all_gene_diff,by = 'Gene_ID',all.y   = T)
expression_diff<-expression_diff[,-c(17:18)]
library(pheatmap)
group<-read.delim('group.txt',row.names = 1)
col<-colorRampPalette(c("blue", "white", "red"))(256)
row.names(expression_diff)<-expression_diff$Gene_ID
expression_diff<-expression_diff[,-1]
df_1<-as.matrix(scale(expression_diff))
diff_heatmap<-pheatmap(df_1,scale = 'row',color=col,show_rownames = F,cellwidth = 25,
                       annotation_col = group,treeheight_row = 30,cluster_cols = F)
diff_heatmap

##GO富集
library(clusterProfiler)
library(org.Mm.eg.db)
GO_up_gene_all<- enrichGO(gene = up_gene_list$Gene_ID,OrgDb= org.Mm.eg.db,
                              keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_up_gene_detail<-as.data.frame(GO_up_gene_all)
row.names(down_gene_list)<-down_gene_list$Gene_ID
GO_down_gene_all<- enrichGO(gene = down_gene_list$Gene_ID,OrgDb= org.Mm.eg.db,
                          keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_down_gene_detail<-as.data.frame(GO_down_gene_all)
#GO富集气泡图(BP)
GO_up_gene_BP_detail<-GO_up_gene_detail[which(GO_up_gene_detail$ONTOLOGY == 'BP'),]
GO_down_gene_BP_detail<-GO_down_gene_detail[which(GO_up_gene_detail$ONTOLOGY == 'BP'),]
GO_up_gene_BP_detail['total']<-155
GO_down_gene_BP_detail['total']<-255
GO_up_gene_BP_detail['Richfactor']<-(GO_up_gene_BP_detail$Count)/(GO_up_gene_BP_detail$total)
GO_down_gene_BP_detail['Richfactor']<-(GO_down_gene_BP_detail$Count)/(GO_down_gene_BP_detail$total)
GO_up_gene_BP_detail<-GO_up_gene_BP_detail[order(GO_up_gene_BP_detail$Richfactor,decreasing = T),]
GO_down_gene_BP_detail<-GO_down_gene_BP_detail[order(GO_down_gene_BP_detail$Richfactor,decreasing = T),]
GO_up_gene_BP_detail_top20<-GO_up_gene_BP_detail[c(1:20),]
GO_down_gene_BP_detail_top20<-GO_down_gene_BP_detail[c(1:20),]
write.csv(GO_up_gene_BP_detail_top20,file = 'Go_up_top20.csv',row.names = GO_up_gene_BP_detail_top20$geneID)
library(ggplot2)
GO_up_gene_BP_detail_top20plot<-ggplot(GO_up_gene_BP_detail_top20,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Biological Process Top20_Up")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_down_gene_BP_detail_top20plot<-ggplot(GO_down_gene_BP_detail_top20,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Biological Process Top20_Down")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_up_gene_BP_detail_top20plot
GO_down_gene_BP_detail_top20plot
ggsave("GO_up_gene_BP_detail_top20plot.png",GO_up_gene_BP_detail_top20plot,width=10,height=8)
ggsave("GO_down_gene_BP_detail_top20plot.png",GO_down_gene_BP_detail_top20plot,width=8,height=8)
##GO富集气泡图(MF)
GO_up_gene_MF_detail<-GO_up_gene_detail[which(GO_up_gene_detail$ONTOLOGY == 'MF'),]
GO_down_gene_MF_detail<-GO_down_gene_detail[which(GO_up_gene_detail$ONTOLOGY == 'MF'),]
GO_up_gene_MF_detail['total']<-155
GO_down_gene_MF_detail['total']<-253
GO_up_gene_MF_detail['Richfactor']<-(GO_up_gene_MF_detail$Count)/(GO_up_gene_MF_detail$total)
GO_down_gene_MF_detail['Richfactor']<-(GO_down_gene_MF_detail$Count)/(GO_down_gene_MF_detail$total)
GO_up_gene_MF_detail<-GO_up_gene_MF_detail[order(GO_up_gene_MF_detail$Richfactor,decreasing = T),]
GO_down_gene_MF_detail<-GO_down_gene_MF_detail[order(GO_down_gene_MF_detail$Richfactor,decreasing = T),]
GO_up_gene_MF_detail_top20<-GO_up_gene_MF_detail[c(1:20),]
GO_down_gene_MF_detail_top20<-GO_down_gene_MF_detail[c(1:20),]
write.csv(GO_up_gene_BP_detail_top20,file = 'Go_up_top20.csv',row.names = GO_up_gene_BP_detail_top20$geneID)
library(ggplot2)
GO_up_gene_MF_detail_top20plot<-ggplot(GO_up_gene_MF_detail_top20,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Molecular Function Top20_Up")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_down_gene_MF_detail_top20plot<-ggplot(GO_down_gene_BP_detail_top20,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Molecular Function Top20_Down")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_up_gene_MF_detail_top20plot
GO_down_gene_MF_detail_top20plot
ggsave("GO_up_gene_MF_detail_top20plot.png",GO_up_gene_MF_detail_top20plot,width=10,height=8)
ggsave("GO_down_gene_MF_detail_top20plot.png",GO_down_gene_MF_detail_top20plot,width=8,height=8)

#KEGG通路富集
library(clusterProfiler)
library(org.Mm.eg.db)
#EMBL编号转换
KEGG_up_geneid<- bitr(unique(up_gene_list$Gene_ID), fromType = "ENSEMBL",
        toType = c( "ENTREZID"),
        OrgDb = org.Mm.eg.db)
KEGG_down_geneid<- bitr(unique(down_gene_list$Gene_ID), fromType = "ENSEMBL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Mm.eg.db)
KEGG_up_gene_all<- enrichKEGG(gene = KEGG_up_geneid$ENTREZID,organism = "mouse",
                              keyType = "kegg",pvalueCutoff = 1,qvalueCutoff = 1,
                              pAdjustMethod = "BH")
KEGG_down_gene_all<- enrichKEGG(gene = KEGG_down_geneid$ENTREZID,organism = "mouse",
                              keyType = "kegg",pvalueCutoff = 1,qvalueCutoff = 1,
                              pAdjustMethod = "BH")
KEGG_up_gene_detail<-as.data.frame(KEGG_up_gene_all)
KEGG_down_gene_detail<-as.data.frame(KEGG_down_gene_all)
browseKEGG(KEGG_up_gene_all,'mmu00830')
##
##分别GO富集
library(clusterProfiler)
library(org.Mm.eg.db)
GO_D10NTC_D10TOA_DEG_up<- enrichGO(gene = D10NTC_D10TOA_DEG_up$id,OrgDb= org.Mm.eg.db,
                          keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D10TOA_DEG_down<- enrichGO(gene = D10NTC_D10TOA_DEG_down$id,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D20TOA_DEG_up<- enrichGO(gene = D10NTC_D20TOA_DEG_up$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D20TOA_DEG_down<- enrichGO(gene = D10NTC_D20TOA_DEG_down$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D30TOA_DEG_up<- enrichGO(gene = D10NTC_D30TOA_DEG_up$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D30TOA_DEG_down<- enrichGO(gene = D10NTC_D30TOA_DEG_down$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D35TOA_DEG_up<- enrichGO(gene = D10NTC_D35TOA_DEG_up$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D35TOA_DEG_down<- enrichGO(gene = D10NTC_D35TOA_DEG_down$Gene_ID,OrgDb= org.Mm.eg.db,
                                   keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
GO_D10NTC_D10TOA_DEG_up_detail<-as.data.frame(GO_D10NTC_D10TOA_DEG_up)
GO_D10NTC_D10TOA_DEG_down_detail<-as.data.frame(GO_D10NTC_D10TOA_DEG_down)
GO_D10NTC_D20TOA_DEG_up_detail<-as.data.frame(GO_D10NTC_D20TOA_DEG_up)
GO_D10NTC_D20TOA_DEG_down_detail<-as.data.frame(GO_D10NTC_D20TOA_DEG_down)
GO_D10NTC_D30TOA_DEG_up_detail<-as.data.frame(GO_D10NTC_D30TOA_DEG_up)
GO_D10NTC_D30TOA_DEG_down_detail<-as.data.frame(GO_D10NTC_D30TOA_DEG_down)
GO_D10NTC_D35TOA_DEG_up_detail<-as.data.frame(GO_D10NTC_D35TOA_DEG_up)
GO_D10NTC_D35TOA_DEG_down_detail<-as.data.frame(GO_D10NTC_D35TOA_DEG_down)
GO_D10NTC_D10TOA_DEG_up_detail['total']<-ifelse(GO_D10NTC_D10TOA_DEG_up_detail$ONTOLOGY == 'BP',60, 59)
GO_D10NTC_D20TOA_DEG_up_detail['total']<-ifelse(GO_D10NTC_D20TOA_DEG_up_detail$ONTOLOGY == 'MF',48,47)
GO_D10NTC_D20TOA_DEG_down_detail['total']<-ifelse(GO_D10NTC_D20TOA_DEG_down_detail$ONTOLOGY == 'BP',129
                                                  ,ifelse(GO_D10NTC_D20TOA_DEG_down_detail == 'CC',132,131))
GO_D10NTC_D30TOA_DEG_up_detail['total']<-ifelse(GO_D10NTC_D30TOA_DEG_up_detail$ONTOLOGY == 'BP',53, 52)
GO_D10NTC_D30TOA_DEG_down_detail['total']<-ifelse(GO_D10NTC_D30TOA_DEG_down_detail$ONTOLOGY == 'BP',101
                                                  ,ifelse(GO_D10NTC_D30TOA_DEG_down_detail == 'CC',98,99))
GO_D10NTC_D35TOA_DEG_up_detail['total']<-ifelse(GO_D10NTC_D35TOA_DEG_up_detail$ONTOLOGY == 'BP',29, 30)
GO_D10NTC_D35TOA_DEG_down_detail['total']<-ifelse(GO_D10NTC_D35TOA_DEG_down_detail$ONTOLOGY == 'BP',49,36)
GO_D10NTC_D10TOA_DEG_up_detail['Richfactor']<-GO_D10NTC_D10TOA_DEG_up_detail$Count/GO_D10NTC_D10TOA_DEG_up_detail$total
GO_D10NTC_D20TOA_DEG_up_detail['Richfactor']<-GO_D10NTC_D20TOA_DEG_up_detail$Count/GO_D10NTC_D20TOA_DEG_up_detail$total
GO_D10NTC_D20TOA_DEG_down_detail['Richfactor']<-GO_D10NTC_D20TOA_DEG_down_detail$Count/GO_D10NTC_D20TOA_DEG_down_detail$total
GO_D10NTC_D30TOA_DEG_up_detail['Richfactor']<-GO_D10NTC_D30TOA_DEG_up_detail$Count/GO_D10NTC_D30TOA_DEG_up_detail$total
GO_D10NTC_D30TOA_DEG_down_detail['Richfactor']<-GO_D10NTC_D30TOA_DEG_down_detail$Count/GO_D10NTC_D35TOA_DEG_down_detail$total
GO_D10NTC_D35TOA_DEG_up_detail['Richfactor']<-GO_D10NTC_D35TOA_DEG_up_detail$Count/GO_D10NTC_D35TOA_DEG_up_detail$total
GO_D10NTC_D35TOA_DEG_down_detail['Richfactor']<-GO_D10NTC_D35TOA_DEG_down_detail$Count/GO_D10NTC_D35TOA_DEG_down_detail$total
GO_D10NTC_D10TOA_DEG_up_detail<-GO_D10NTC_D10TOA_DEG_up_detail[order(GO_D10NTC_D10TOA_DEG_up_detail$Richfactor,decreasing = T),]
GO_D10NTC_D20TOA_DEG_up_detail<-GO_D10NTC_D20TOA_DEG_up_detail[order(GO_D10NTC_D20TOA_DEG_up_detail$Richfactor,decreasing = T),]
GO_D10NTC_D20TOA_DEG_down_detail<-GO_D10NTC_D20TOA_DEG_down_detail[order(GO_D10NTC_D20TOA_DEG_down_detail$Richfactor,decreasing = T),]
GO_D10NTC_D30TOA_DEG_up_detail<-GO_D10NTC_D30TOA_DEG_up_detail[order(GO_D10NTC_D30TOA_DEG_up_detail$Richfactor,decreasing = T),]
GO_D10NTC_D30TOA_DEG_down_detail<-GO_D10NTC_D30TOA_DEG_down_detail[order(GO_D10NTC_D30TOA_DEG_down_detail$Richfactor,decreasing = T),]
GO_D10NTC_D35TOA_DEG_up_detail<-GO_D10NTC_D35TOA_DEG_up_detail[order(GO_D10NTC_D35TOA_DEG_up_detail$Richfactor,decreasing = T),]
GO_D10NTC_D35TOA_DEG_down_detail<-GO_D10NTC_D35TOA_DEG_up_detail[order(GO_D10NTC_D35TOA_DEG_down_detail$Richfactor,decreasing = T),]
write.csv(file = 'GO_D10NTC_D10TOA_DEG_up_detail.csv',GO_D10NTC_D10TOA_DEG_up_detail)
write.csv(file = 'GO_D10NTC_D20TOA_DEG_up_detail.csv',GO_D10NTC_D20TOA_DEG_up_detail)
write.csv(file = 'GO_D10NTC_D30TOA_DEG_up_detail.csv',GO_D10NTC_D30TOA_DEG_up_detail)
write.csv(file = 'GO_D10NTC_D35TOA_DEG_up_detail.csv',GO_D10NTC_D35TOA_DEG_up_detail)
#GO柱状图汇总
library(ggplot2)
GO_D10NTC_D10TOA_DEG_up_detail_plot<-ggplot(GO_D10NTC_D10TOA_DEG_up_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=4,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D10TOA_DEG_up_detail_plot
GO_D10NTC_D20TOA_DEG_up_detail_plot<-ggplot(GO_D10NTC_D20TOA_DEG_up_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=5,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D20TOA_DEG_up_detail_plot
GO_D10NTC_D20TOA_DEG_down_detail_plot<-ggplot(GO_D10NTC_D20TOA_DEG_down_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=3,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D20TOA_DEG_down_detail_plot
GO_D10NTC_D30TOA_DEG_up_detail_plot<-ggplot(GO_D10NTC_D30TOA_DEG_up_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=5,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D30TOA_DEG_up_detail_plot
GO_D10NTC_D30TOA_DEG_down_detail_plot<-ggplot(GO_D10NTC_D30TOA_DEG_down_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=3,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D30TOA_DEG_down_detail_plot
GO_D10NTC_D35TOA_DEG_up_detail_plot<-ggplot(GO_D10NTC_D35TOA_DEG_up_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=5,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D35TOA_DEG_up_detail_plot
GO_D10NTC_D35TOA_DEG_down_detail_plot<-ggplot(GO_D10NTC_D35TOA_DEG_down_detail,aes(x=ID,y=Count,fill=pvalue))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+geom_text(aes(label=Count),size=5,vjust= -0.5)+
  ylab('Gene Count')+xlab('GO ID')+theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size= 12),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=8,colour = 'black',angle = 90),legend.text = element_text(size=13,colour = 'black'))
GO_D10NTC_D35TOA_DEG_down_detail_plot
ggsave('GO_D10NTC_D10TOA_up.tiff',GO_D10NTC_D10TOA_DEG_up_detail_plot,height = 6,width = 8)
ggsave('GO_D10NTC_D20TOA_up.tiff',GO_D10NTC_D20TOA_DEG_up_detail_plot,height = 6,width = 8)
ggsave('GO_D10NTC_D20TOA_down.tiff',GO_D10NTC_D20TOA_DEG_down_detail_plot,height = 6,width = 12)
ggsave('GO_D10NTC_D30TOA_up.tiff',GO_D10NTC_D30TOA_DEG_up_detail_plot,height = 6,width = 8)
ggsave('GO_D10NTC_D30TOA_down.tiff',GO_D10NTC_D30TOA_DEG_down_detail_plot,height = 6,width = 12)
ggsave('GO_D10NTC_D35TOA_up.tiff',GO_D10NTC_D35TOA_DEG_up_detail_plot,height = 6,width = 8)
ggsave('GO_D10NTC_D35TOA_down.tiff',GO_D10NTC_D35TOA_DEG_down_detail_plot,height = 6,width = 8)
#GO气泡图汇总
library(ggplot2)
GO_D10NTC_D10TOA_DEG_up_top20plot<-ggplot(GO_D10NTC_D10TOA_DEG_up_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D10TOA_DEG_up_top20plot
GO_D10NTC_D20TOA_DEG_up_top20plot<-ggplot(GO_D10NTC_D20TOA_DEG_up_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D20TOA_DEG_up_top20plot
GO_D10NTC_D20TOA_DEG_down_top20plot<-ggplot(GO_D10NTC_D20TOA_DEG_down_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D20TOA_DEG_down_top20plot
GO_D10NTC_D30TOA_DEG_up_top20plot<-ggplot(GO_D10NTC_D30TOA_DEG_up_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D30TOA_DEG_up_top20plot
GO_D10NTC_D30TOA_DEG_down_top20plot<-ggplot(GO_D10NTC_D30TOA_DEG_down_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D30TOA_DEG_down_top20plot
GO_D10NTC_D35TOA_DEG_up_top20plot<-ggplot(GO_D10NTC_D35TOA_DEG_up_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D35TOA_DEG_up_top20plot
GO_D10NTC_D35TOA_DEG_down_top20plot<-ggplot(GO_D10NTC_D35TOA_DEG_down_detail[c(1:20),],aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_D10NTC_D35TOA_DEG_down_top20plot
ggsave('GO_D10NTC_D10TOA_DEG_up_top20plot.tiff',GO_D10NTC_D10TOA_DEG_up_top20plot,width = 12,height = 8)
ggsave('GO_D10NTC_D20TOA_DEG_up_top20plot.tiff',GO_D10NTC_D20TOA_DEG_up_top20plot,width = 10,height = 8)
ggsave('GO_D10NTC_D20TOA_DEG_down_top20plot.tiff',GO_D10NTC_D20TOA_DEG_down_top20plot,width = 8,height = 8)
ggsave('GO_D10NTC_D30TOA_DEG_up_top20plot.tiff',GO_D10NTC_D30TOA_DEG_up_top20plot,width = 10,height = 8)
ggsave('GO_D10NTC_D30TOA_DEG_down_top20plot.tiff',GO_D10NTC_D30TOA_DEG_down_top20plot,width = 8,height = 8)
ggsave('GO_D10NTC_D35TOA_DEG_up_top20plot.tiff',GO_D10NTC_D35TOA_DEG_up_top20plot,width = 10,height = 8)
ggsave('GO_D10NTC_D35TOA_DEG_down_top20plot.tiff',GO_D10NTC_D35TOA_DEG_down_top20plot,width = 10,height = 8)
##趋势分析
##分析anti表中数据
anti<-read.delim('anti.txt')
library(clusterProfiler)
library(org.Mm.eg.db)
anti<-bitr(unique(anti$Gene_ID), fromType = "ENSEMBL",toType = c( "SYMBOL"),OrgDb = org.Mm.eg.db)
anti['Gene_ID']<-anti$ENSEMBL
anti_1<-merge(anti,expression_test_fpkm,by = 'Gene_ID',all.x = T)
anti_1<-anti_1[,-2]
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
anti_2<-anti_1[,-c(1,3:22)]
anti_2<-anti_2[order(anti_2$D10TOA,decreasing = T),]
library(reshape2)
anti_2.1<-melt(anti_2,variable.name = 'Group',value.name = 'rate')
anti_2.1['label']<-ifelse(anti_2.1$Group == 'D35TOA',anti_2$SYMBOL,'')
library(ggplot2)
anti_2.1_plot<-ggplot(anti_2.1,aes(x=Group,y=rate,colour=SYMBOL,group=SYMBOL))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide gene plot",x='Group',y='log2(V1/V0)')+geom_text(aes(label=label),size=3,vjust=1)+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text('Antimicrobial peptide gene'),
        legend.text =element_text(size=11,colour = 'black') )
anti_2.1_plot
ggsave('anti_2plot.png',anti_2.1_plot,width = 8,height = 8)
#正趋势基因折线图
anti2.2<-anti_2
anti2.2['stat']<-ifelse(anti2.2$D10TOA >0, 'up','NA')
p<-grep('NA',anti2.2$stat)
anti2.2<-anti2.2[-p,]
library(reshape2)
anti2.2<-melt(anti2.2,variable.name = 'Group',value.name = 'rate')
anti2.2['label']<-ifelse(anti2.2$Group == 'D35TOA',anti2.2$SYMBOL, '')
library(ggplot2)
anti2.2_plot<-ggplot(anti2.2,aes(x=Group,y=rate,colour=SYMBOL,group=SYMBOL))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide gene plot",x='Group',y=expression(log[2](Fold-Change)))+geom_text(aes(label = label),size=3,vjust=1)+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text('Antimicrobial peptide gene'),
        legend.text =element_text(size=11,colour = 'black') )
anti2.2_plot
ggsave('anti_2.2plot.png',anti2.2_plot,width = 8,height = 8)

anti_3<-anti_2
row.names(anti_3)<-anti_3$SYMBOL
anti_3<-anti_3[,-1]
anti_3.1<-as.matrix(anti_3)
library(pheatmap)
anti_heatmap<-pheatmap(anti_3,cutree_rows = 6,cluster_cols = F)
anti_heatmap

##
##暂时完成
##
##
##又开始了!开不开心？？？
##以D10PBS为空白进行重分析
##
##DEseq2差异基因分组分析
#分别整理各组分析表
expression_all<-read.delim('Expression_all.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                           check.names = F)
##分析D10PBS与饲喂d10，d20，d30，d35差异基因情况
#摘取不同组counts数据
D10PBS_D20TOA_count<-expression_all[,c(8,10,12,32,34,36)]
D10PBS_D30TOA_count<-expression_all[,c(8,10,12,50,52,54)]
D10PBS_D35TOA_count<-expression_all[,c(8,10,12,56,58,60)]
#编写样本组
condition_4 <- factor(c(rep("D10PBS",3),rep("D20TOA",3)), levels = c("D10PBS","D20TOA"))#一定为factor格式
condition_5<- factor(c(rep("D10PBS",3),rep("D30TOA",3)), levels = c("D10PBS","D30TOA"))
condition_6 <- factor(c(rep("D10PBS",3),rep("D35TOA",3)), levels = c("D10PBS","D35TOA"))
#结合形成样本组数据框
coldata_4<-data.frame(row.names = colnames(D10PBS_D20TOA_count),condition_4)
coldata_5<-data.frame(row.names = colnames(D10PBS_D30TOA_count),condition_5)
coldata_6<-data.frame(row.names = colnames(D10PBS_D35TOA_count),condition_6)
#dds=DESeqDataSet Object构建dss模型
library(DESeq2)
#D10PBS与D20TOA
dds_4<- DESeqDataSetFromMatrix(D10PBS_D20TOA_count, coldata_4, design= ~ condition_4)
dds_4 <- DESeq(dds_4)
#总体结果查看
res_4<- results(dds_4, contrast=c("condition_4", "D20TOA", "D10PBS"))
res_4<-na.omit(res_4)
res_4_detail<-as.data.frame(res_4)
write.csv(res_4,file="D10PBS_D20TOA_results.csv")
#差异基因输出
D10PBS_D20TOA_DEG_up<-subset(res_4, pvalue < 0.05 & log2FoldChange >1)#DEG上调
D10PBS_D20TOA_DEG_down<-subset(res_4, pvalue < 0.05 & log2FoldChange < -1)#DEG下调
#D10PBS与D30TOA
dds_5<- DESeqDataSetFromMatrix(D10PBS_D30TOA_count, coldata_5, design= ~ condition_5)
dds_5 <- DESeq(dds_5)
#总体结果查看
res_5<- results(dds_5, contrast=c("condition_5", "D30TOA", "D10PBS"))
res_5<-na.omit(res_5)
#差异基因输出
D10PBS_D30TOA_DEG_up<-subset(res_5, pvalue < 0.05 & log2FoldChange > 1)#DEG上调
D10PBS_D30TOA_DEG_down<-subset(res_5, pvalue < 0.05 & log2FoldChange < -1)#DEG下调
#D10PBS与D35TOA
dds_6<- DESeqDataSetFromMatrix(D10PBS_D35TOA_count, coldata_6, design= ~ condition_6)
dds_6<- DESeq(dds_6)
#总体结果查看
res_6<-results(dds_6, contrast=c("condition_6", "D35TOA", "D10PBS"))
res_6<-na.omit(res_6)
#差异基因输出
D10PBS_D35TOA_DEG_up<-subset(res_6, pvalue < 0.05 & log2FoldChange >1)#DEG上调
D10PBS_D35TOA_DEG_down<-subset(res_6, pvalue < 0.05 & log2FoldChange< -1)#DEG下调
#火山图
library(readxl)
D10PBS_D10TOA_DEG<-read_xlsx('D10PBS_vs_D10TOA.DESeq.xlsx')
D10PBS_D20TOA_DEG<-as.data.frame(res_4)
D10PBS_D30TOA_DEG<-as.data.frame(res_5)
D10PBS_D35TOA_DEG<-as.data.frame(res_6)
#添加上升下降趋势
D10PBS_D10TOA_DEG<-D10PBS_D10TOA_DEG[,c(1,3:7)]
D10PBS_D10TOA_DEG<-D10PBS_D10TOA_DEG[-(which(D10PBS_D10TOA_DEG$baseMeanA == 0 | D10PBS_D10TOA_DEG$baseMeanB ==0)),]#删除带INF值
D10PBS_D10TOA_DEG$log2FoldChange<-as.numeric(D10PBS_D10TOA_DEG$log2FoldChange)
D10PBS_D10TOA_DEG['significant']<-ifelse(
  D10PBS_D10TOA_DEG$pval < 0.05 & D10PBS_D10TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10PBS_D10TOA_DEG$pval < 0.05 & D10PBS_D10TOA_DEG$log2FoldChange < -1,"down","no"))
D10PBS_D20TOA_DEG['significant']<-ifelse(
  D10PBS_D20TOA_DEG$pvalue < 0.05 & D10PBS_D20TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10PBS_D20TOA_DEG$pvalue < 0.05 & D10PBS_D20TOA_DEG$log2FoldChange < -1,"down","no"))
D10PBS_D30TOA_DEG['significant']<-ifelse(
  D10PBS_D30TOA_DEG$pvalue < 0.05 & D10PBS_D30TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10PBS_D30TOA_DEG$pvalue < 0.05 & D10PBS_D30TOA_DEG$log2FoldChange < -1,"down","no"))
D10PBS_D35TOA_DEG['significant']<-ifelse(
  D10PBS_D35TOA_DEG$pvalue < 0.05 & D10PBS_D35TOA_DEG$log2FoldChange >1,"up",
  ifelse(D10PBS_D35TOA_DEG$pvalue < 0.05 & D10PBS_D35TOA_DEG$log2FoldChange < -1,"down","no"))
#ggplot2火山图可视化
library(ggplot2)
D10PBSvsD10TOA_plot<-ggplot(D10PBS_D10TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pval)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10NTC_vs_D10TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),
        legend.text = element_text(size=15,colour = 'black'),legend.title = element_text(size=15,colour = 'black'))
D10PBSvsD20TOA_plot<-ggplot(D10PBS_D20TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10PBS_vs_D20TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10PBSvsD30TOA_plot<-ggplot(D10PBS_D30TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10PBS_vs_D30TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10PBSvsD35TOA_plot<-ggplot(D10PBS_D35TOA_DEG,aes(x=log2FoldChange,y=-1*log10(pvalue)))+geom_point(aes(color=significant),size=2)+
  xlim(-5,5)+ylim(0,7)+labs(title="D10PBS_vs_D35TOA Volcano Plot",
                            x=expression(log[2](FC)),y=expression(-log[10](p-value)))+
  scale_color_manual(values =c("blue",'gray','red'))+geom_hline(yintercept=1.3,linetype=2,col="black")+
  geom_vline(xintercept=c(-1,1),linetype=2,col="black")+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title = element_text(size=15,colour = 'black'))
D10PBSvsD10TOA_plot
D10PBSvsD20TOA_plot
D10PBSvsD30TOA_plot
D10PBSvsD35TOA_plot
ggsave("D10PBSvsD10TOA_volcano plot.png",D10PBSvsD10TOA_plot,width=8,height=8)
ggsave("D10PBSvsD20TOA_volcano plot.png",D10PBSvsD20TOA_plot,width=8,height=8)
ggsave("D10PBSvsD30TOA_volcano plot.png",D10PBSvsD30TOA_plot,width=8,height=8)
ggsave("D10PBSvsD35TOA_volcano plot.png",D10PBSvsD35TOA_plot,width=8,height=8)

#得到的上升下降的DEG转换为数据框
D10PBS_D20TOA_DEG_up<-as.data.frame(D10PBS_D20TOA_DEG_up)
D10PBS_D20TOA_DEG_down<-as.data.frame(D10PBS_D20TOA_DEG_down)
D10PBS_D30TOA_DEG_up<-as.data.frame(D10PBS_D30TOA_DEG_up)
D10PBS_D30TOA_DEG_down<-as.data.frame(D10PBS_D30TOA_DEG_down)
D10PBS_D35TOA_DEG_up<-as.data.frame(D10PBS_D35TOA_DEG_up)
D10PBS_D35TOA_DEG_down<-as.data.frame(D10PBS_D35TOA_DEG_down)
D10PBS_D20TOA_DEG_up['Gene_ID']<-row.names(D10PBS_D20TOA_DEG_up)
D10PBS_D20TOA_DEG_down['Gene_ID']<-row.names(D10PBS_D20TOA_DEG_down)
D10PBS_D30TOA_DEG_up['Gene_ID']<-row.names(D10PBS_D30TOA_DEG_up)
D10PBS_D30TOA_DEG_down['Gene_ID']<-row.names(D10PBS_D30TOA_DEG_down)
D10PBS_D35TOA_DEG_up['Gene_ID']<-row.names(D10PBS_D35TOA_DEG_up)
D10PBS_D35TOA_DEG_down['Gene_ID']<-row.names(D10PBS_D35TOA_DEG_down)
#获得symbol列表
expression_all['Gene_ID']<-row.names(expression_all)
gene_detail<-expression_all[,c(67:75,77)]
#merge函数结合数据框
D10PBS_D20TOA_DEG_up<-merge(D10PBS_D20TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10PBS_D20TOA_DEG_down<-merge(D10PBS_D20TOA_DEG_down,gene_detail,by='Gene_ID',all.x=T)
D10PBS_D30TOA_DEG_up<-merge(D10PBS_D30TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10PBS_D30TOA_DEG_down<-merge(D10PBS_D30TOA_DEG_down,gene_detail,by='Gene_ID',all.x=T)
D10PBS_D35TOA_DEG_up<-merge(D10PBS_D35TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
D10PBS_D35TOA_DEG_down<-merge(D10PBS_D35TOA_DEG_up,gene_detail,by='Gene_ID',all.x=T)
#输出数据csv文件
write.csv(D10PBS_D20TOA_DEG_up,file = 'D10PBS_D20TOA_DEG_up_gene.csv',
          row.names = D10PBS_D20TOA_DEG_up$Gene_ID)
write.csv(D10PBS_D20TOA_DEG_down,file = 'D10PBS_D20TOA_DEG_down_gene.csv',
          row.names = D10PBS_D20TOA_DEG_down$Gene_ID)
write.csv(D10PBS_D30TOA_DEG_up,file = 'D10PBS_D30TOA_DEG_up_gene.csv',
          row.names = D10PBS_D30TOA_DEG_up$Gene_ID)
write.csv(D10PBS_D30TOA_DEG_down,file = 'D10PBS_D30TOA_DEG_down_gene.csv',
          row.names = D10PBS_D30TOA_DEG_down$Gene_ID)
write.csv(D10PBS_D35TOA_DEG_up,file = 'D10PBS_D35TOA_DEG_up_gene.csv',
          row.names = D10PBS_D20TOA_DEG_up$Gene_ID)
write.csv(D10PBS_D35TOA_DEG_down,file = 'D10PBS_D35TOA_DEG_down_gene.csv',
          row.names = D10PBS_D35TOA_DEG_down$Gene_ID)
#柱状图统计
detail_PBS<-read.delim('detail_PBS.txt')
library(reshape2)
detail_2<-melt(detail_PBS,variable.name = 'Type',value.name = 'Count')
library(ggplot2)
diff.stat.plot<-ggplot(detail_2,aes(x=Group,y=Count,fill=Type))+
  geom_bar(stat ="identity",width = 0.8,position = "dodge")+coord_flip()+ 
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=12),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text = element_text(size=15,colour = 'black'),
        legend.title =element_text(size = 15,color = 'black') )
diff.stat.plot
ggsave('diff_stat.png',diff.stat.plot,width = 8,height = 8)
#查找4组相同上升，下降基因(UpSetR)
library(UpSetR)
library(readxl)
D10PBS_D10TOA_DEG_up<-read_xlsx('D10PBS_vs_D10TOA.DESeq.Up.xlsx')
D10PBS_D10TOA_DEG_down<-read_xlsx('D10PBS_vs_D10TOA.DESeq.Down.xlsx')
listinput_up_pbs<-list(D10PBS_D10TOA_DEG_up= D10PBS_D10TOA_DEG_up$id,
                   D10PBS_D20TOA_DEG_up = D10PBS_D20TOA_DEG_up$Gene_ID,
                   D10PBS_D30TOA_DEG_up = D10PBS_D30TOA_DEG_up$Gene_ID,
                   D10PBS_D35TOA_DEG_up = D10PBS_D35TOA_DEG_up$Gene_ID)
listinput_down_pbs<-list(D10PBS_D10TOA_DEG_down =D10PBS_D10TOA_DEG_down$id,
                     D10PBS_D20TOA_DEG_down =D10PBS_D20TOA_DEG_down$Gene_ID,
                     D10PBS_D30TOA_DEG_down =D10PBS_D30TOA_DEG_down$Gene_ID,
                     D10PBS_D35TOA_DEG_down =D10PBS_D35TOA_DEG_down$Gene_ID)
listinput_up_plot_pbs<- upset(fromList(listinput_up_pbs),nsets = 5, order.by = "freq")
listinput_down_plot_pbs<-upset(fromList(listinput_down_pbs),nsets = 5, order.by = "freq")
listinput_up_plot_pbs
listinput_down_plot_pbs
##热图
#获得所有组上调基因EMBL编号
D10PBS_D10TOA_DEG_up['Gene_ID']<-D10PBS_D10TOA_DEG_up$id
D10PBS_D10TOA_DEG_down['Gene_ID']<-D10PBS_D10TOA_DEG_down$id
up_gene_list_PBS<-merge(D10PBS_D20TOA_DEG_up,D10PBS_D10TOA_DEG_up,by = 'Gene_ID',all = T )
up_gene_list_PBS<-merge(up_gene_list_PBS,D10PBS_D30TOA_DEG_up,by = 'Gene_ID',all = T)
up_gene_list_PBS<-merge(up_gene_list_PBS,D10PBS_D35TOA_DEG_up,by = 'Gene_ID',all = T)
row.names(up_gene_list_PBS)<-up_gene_list_PBS$Gene_ID
up_gene_list_PBS<-up_gene_list_PBS[,1]
up_gene_list_PBS<-as.data.frame(up_gene_list_PBS)
up_gene_list_PBS['Gene_ID']<-up_gene_list_PBS$up_gene_list_PBS
#获得所有组下调基因EMBL编号
down_gene_list_PBS<-merge(D10PBS_D20TOA_DEG_down,D10PBS_D10TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list_PBS<-merge(down_gene_list_PBS,D10PBS_D30TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list_PBS<-merge(down_gene_list_PBS,D10PBS_D35TOA_DEG_down,by = 'Gene_ID',all = T)
down_gene_list_PBS<-down_gene_list_PBS[,1]
down_gene_list_PBS<-as.data.frame(down_gene_list_PBS)
down_gene_list_PBS['Gene_ID']<-down_gene_list_PBS$down_gene_list_PBS
all_gene_diff_PBS<-merge(up_gene_list_PBS,down_gene_list_PBS,by = 'Gene_ID',all = T)
#获得所有基因fpkm值
expression_test_PBSfpkm<-expression_all[,c(9,11,13,15,17,19,33,35,37,51,53,55,57,59,61,77)]
#get diff gene fpkm value
expression_diff_PBS<-merge(expression_test_PBSfpkm,all_gene_diff_PBS,by = 'Gene_ID',all.y   = T)
expression_diff_PBS<-expression_diff_PBS[,-c(17:18)]
library(pheatmap)
group<-read.delim('group.txt',row.names = 1)
col<-colorRampPalette(c("blue", "white", "red"))(256)
row.names(expression_diff_PBS)<-expression_diff_PBS$Gene_ID
expression_diff_PBS<-expression_diff_PBS[,-1]
df_2<-as.matrix(scale(expression_diff_PBS))
diff_heatmap_PBS<-pheatmap(df_2,scale = 'row',color=col,show_rownames = F,cellwidth = 25,
                       annotation_col = group,treeheight_row = 30,cluster_cols = F)
diff_heatmap_PBS
##GO富集
library(clusterProfiler)
library(org.Mm.eg.db)
#EMBL编号转换
row.names(up_gene_list_PBS)<-up_gene_list_PBS$Gene_ID
up_gene_list_PBS['ensembl_id']<-up_gene_list_PBS$Gene_ID
Entrez_gene_id<-toTable(org.Mm.egENSEMBL)
up_gene_list_GO_PBS<-merge(Entrez_gene_id,up_gene_list_PBS,by='ensembl_id',all.x = T)
up_gene_list_GO_PBS<-na.omit(up_gene_list_GO_PBS)
row.names(up_gene_list_PBS)<-up_gene_list_PBS$gene_id
GO_up_gene_all_PBS<- enrichGO(gene = up_gene_list_PBS$Gene_ID,OrgDb= org.Mm.eg.db,
                          keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.1)
GO_up_gene_detail_PBS<-as.data.frame(GO_up_gene_all_PBS)
row.names(down_gene_list_PBS)<-down_gene_list_PBS$Gene_ID
GO_down_gene_all_PBS<- enrichGO(gene = down_gene_list_PBS$Gene_ID,OrgDb= org.Mm.eg.db,
                            keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.1)
GO_down_gene_detail_PBS<-as.data.frame(GO_down_gene_all_PBS)
#GO富集气泡图
GO_up_gene_BP_detail_PBS<-GO_up_gene_detail_PBS[which(GO_up_gene_detail_PBS$ONTOLOGY == 'BP'),]
GO_down_gene_BP_detail_PBS<-GO_down_gene_detail_PBS[which(GO_down_gene_detail_PBS$ONTOLOGY == 'BP'),]
GO_up_gene_BP_detail_PBS['total']<-498
GO_down_gene_BP_detail_PBS['total']<-807
GO_up_gene_BP_detail_PBS['Richfactor']<-(GO_up_gene_BP_detail_PBS$Count)/(GO_up_gene_BP_detail_PBS$total)
GO_down_gene_BP_detail_PBS['Richfactor']<-(GO_down_gene_BP_detail_PBS$Count)/(GO_down_gene_BP_detail_PBS$total)
GO_up_gene_BP_detail_PBS<-GO_up_gene_BP_detail_PBS[order(GO_up_gene_BP_detail_PBS$Richfactor,decreasing = T),]
GO_down_gene_BP_detail_PBS<-GO_down_gene_BP_detail_PBS[order(GO_down_gene_BP_detail_PBS$Richfactor,decreasing = T),]
GO_up_gene_BP_detail_top20_PBS<-GO_up_gene_BP_detail_PBS[c(1:20),]
GO_down_gene_BP_detail_top20_PBS<-GO_down_gene_BP_detail_PBS[c(1:20),]
library(ggplot2)
GO_up_gene_BP_detail_top20_PBS_plot<-ggplot(GO_up_gene_BP_detail_top20_PBS,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Biological Process Top20_Up")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_down_gene_BP_detail_top20_PBS_plot<-ggplot(GO_down_gene_BP_detail_top20_PBS,aes(-1*log10(Richfactor),Description))+
  geom_point(aes(size=Count,color= -1*log10(pvalue)))+scale_color_gradient(low = "blue",high = "red")+labs(title="Biological Process Top20_Down")+
  theme(plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
        axis.text.x = element_text(size=12,colour = 'black'))
GO_up_gene_BP_detail_top20_PBS_plot
GO_down_gene_BP_detail_top20_PBS_plot
ggsave("GO_up_gene_BP_detail_top20_PBS_plot.png",GO_up_gene_BP_detail_top20_PBS_plot,width=15,height=8)
ggsave("GO_down_gene_BP_detail_top20_PBS_plot.png",GO_down_gene_BP_detail_top20_PBS_plot,width=8,height=8)
#KEGG通路富集
library(clusterProfiler)
library(org.Mm.eg.db)
KEGG_all_geneid<-bitr(unique(expression_diff$Gene_ID), fromType = "ENSEMBL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Mm.eg.db)
KEGG_up_geneid_PBS<- bitr(unique(up_gene_list_PBS$Gene_ID), fromType = "ENSEMBL",
                      toType = c( "ENTREZID"),
                      OrgDb = org.Mm.eg.db)
KEGG_down_geneid_PBs<- bitr(unique(down_gene_list_PBS$Gene_ID), fromType = "ENSEMBL",
                        toType = c( "ENTREZID"),
                        OrgDb = org.Mm.eg.db)
KEGG_up_gene_all<- enrichKEGG(gene = KEGG_up_geneid_PBS$ENTREZID,organism = "mouse",
                              keyType = "kegg",pvalueCutoff = 1,qvalueCutoff = 1,
                              pAdjustMethod = "BH")
KEGG_down_gene_all<- enrichKEGG(gene = KEGG_down_geneid_PBs$ENTREZID,organism = "mouse",
                                keyType = "kegg",pvalueCutoff = 1,qvalueCutoff = 1,
                                pAdjustMethod = "BH")
KEGG_up_gene_detail<-as.data.frame(KEGG_up_gene_all)
KEGG_down_gene_detail<-as.data.frame(KEGG_down_gene_all)
browseKEGG(KEGG_up_gene_all,'mmu00830')
##趋势分析
##分析anti表中数据
anti_3<-merge(anti,expression_test_PBSfpkm,by = 'Gene_ID',all.x = T)
anti_3<-anti_3[,-3]
anti_3['D10PBS_mean']<-apply(anti_3[,c(3:5)], 1, mean)
anti_3['D10TOA_mean']<-apply(anti_3[,c(6:8)], 1, mean)
anti_3['D20TOA_mean']<-apply(anti_3[,c(9:11)], 1, mean)
anti_3['D30TOA_mean']<-apply(anti_3[,c(12:14)], 1, mean)
anti_3['D35TOA_mean']<-apply(anti_3[,c(15:17)], 1, mean)
anti_3['D10PBS']<-0
anti_3['D10TOA']<-log2(anti_3$D10TOA_mean/anti_3$D10PBS_mean)
anti_3['D20TOA']<-log2(anti_3$D20TOA_mean/anti_3$D10PBS_mean)
anti_3['D30TOA']<-log2(anti_3$D30TOA_mean/anti_3$D10PBS_mean)
anti_3['D35TOA']<-log2(anti_3$D35TOA_mean/anti_3$D10PBS_mean)
anti_4<-anti_3[,-c(1,3:22)]
anti_4<-anti_4[-c(1,3),]
library(reshape2)
anti_4.1<-melt(anti_4,variable.name = 'Group',value.name = 'rate')
library(ggplot2)
anti_4.1_plot<-ggplot(anti_4.1,aes(x=Group,y=rate,colour=Name,group=Name))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide gene plot",x='Group',y='log2(V1/V0)')+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text('Antimicrobial peptide gene'),
        legend.text =element_text(size=11,colour = 'black') )
anti_4.1_plot
ggsave('anti_4plot.png',anti_4.1_plot,width = 8,height = 8)
#正趋势4基因折线图
anti_4<-anti_4[order(anti_4$D20TOA,decreasing = T),]
anti_4.2<-anti_4[c(1:4),]
library(reshape2)
anti_4.2<-melt(anti_4.2,variable.name = 'Group',value.name = 'rate')
library(ggplot2)
anti_4.2_plot<-ggplot(anti_4.2,aes(x=Group,y=rate,colour=Name,group=Name))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide gene plot",x='Group',y='log2(V1/V0)')+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text('Antimicrobial peptide gene'),
        legend.text =element_text(size=11,colour = 'black') )
anti_4.2_plot
ggsave('anti_4.2plot.png',anti_4.2_plot,width = 8,height = 8)
##整体已运行完

##
library(clusterProfiler)
library(org.Mm.eg.db)
D10NTC_D30TOA_DEG_gene_SYMBOL<-bitr(unique(row.names(D10NTC_D30TOA_DEG)), fromType = "ENSEMBL",
                      toType = c( "SYMBOL"),OrgDb = org.Mm.eg.db)

D10NTC_D30TOA_DEG['ENSEMBL']<-row.names(D10NTC_D30TOA_DEG)
D10NTC_D30TOA_DEG<-D10NTC_D30TOA_DEG[,-8]
D10NTC_D30TOA_DEG_1<-merge(D10NTC_D30TOA_DEG,D10NTC_D30TOA_DEG_gene_SYMBOL,by='ENSEMBL',all.x=T)
D10NTC_D30TOA_DEG_1<-na.omit(D10NTC_D30TOA_DEG_1)
write.csv(file = 'D10NTC_D30TOA_DEG_1.csv',D10NTC_D30TOA_DEG_1,row.names = D10NTC_D30TOA_DEG_1$ ENSEMBL)

##qPCR结果趋势图
qPCR_data<-read.delim('qPCR_data.txt')
qPCR_data['Defa39_gene']<-log2(apply(qPCR_data[,c(2:4)],1,mean))
qPCR_data['Defa35_gene']<-log2(apply(qPCR_data[,c(5:7)],1,mean))
qPCR_data['Defa38_gene']<-log2(apply(qPCR_data[,c(8:10)],1,mean))
qPCR_data['Reg1_gene']<-log2(apply(qPCR_data[,c(11:13)],1,mean))
qPCR_data['Reg3a_gene']<-log2(apply(qPCR_data[,c(14:16)],1,mean))
qPCR_data<-qPCR_data[,c(1,17:21)]
library(reshape2)
qPCR_data_1<-melt(qPCR_data,variable.name = 'Group',value.name = 'Rate')
library(ggplot2)
qPCR_data_plot<-ggplot(qPCR_data_1,aes(x=Name,y=Rate,colour=Group,group=Group))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide qPCR result",x='',y=expression(log[2] (Fold-Change)))+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text =element_text(size=13,colour = 'black'),
        legend.title =element_text(size = 13,color = 'black')  )
qPCR_data_plot
qPCR_data_2<-melt(qPCR_data[,-5],variable.name = 'Group',value.name = 'Rate')
library(ggplot2)
qPCR_data_plot_1<-ggplot(qPCR_data_2,aes(x=Name,y=Rate,colour=Group,group=Group))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(title="Antimicrobial peptide qPCR result",x='',y=expression(log[2] (Fold-Change)))+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.text =element_text(size=13,colour = 'black'),
        legend.title =element_text(size = 13,color = 'black')  )
qPCR_data_plot_1
ggsave('qPCR_data_plot.tiff',qPCR_data_plot,width = 8,height = 8)
ggsave('qPCR_data_plot_1.tiff',qPCR_data_plot_1,width = 8,height = 8)

##测试
tight_protein<-read.delim('tight_protein.txt')
tight_protein<-merge(tight_protein,expression_test_fpkm,by = 'Gene_ID',all.x = T)
tight_protein<-tight_protein[,-2]
tight_protein['D10NTC_mean']<-apply(tight_protein[,c(2:4)], 1, mean)
tight_protein['D10TOA_mean']<-apply(tight_protein[,c(5:7)], 1, mean)
tight_protein['D20TOA_mean']<-apply(tight_protein[,c(8:10)], 1, mean)
tight_protein['D30TOA_mean']<-apply(tight_protein[,c(11:13)], 1, mean)
tight_protein['D35TOA_mean']<-apply(tight_protein[,c(14:16)], 1, mean)
tight_protein['D10NTC']<-0
tight_protein['D10TOA']<-log2(tight_protein$D10TOA_mean/tight_protein$D10NTC_mean)
tight_protein['D20TOA']<-log2(tight_protein$D20TOA_mean/tight_protein$D10NTC_mean)
tight_protein['D30TOA']<-log2(tight_protein$D30TOA_mean/tight_protein$D10NTC_mean)
tight_protein['D35TOA']<-log2(tight_protein$D35TOA_mean/tight_protein$D10NTC_mean)
tight_protein<-tight_protein[,c(1,22:26)]
library(reshape2)
tight_protein_1<-melt(tight_protein,variable.name = 'Group',value.name = 'rate')
library(ggplot2)
tight_plot<-ggplot(tight_protein_1,aes(x=Group,y=rate,colour=Gene_ID,group=Gene_ID))+geom_line()+
  geom_point()+geom_hline(yintercept=0,linetype=2,col="black")+
  labs(x='Group',y='log2(V1/V0)')+
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),plot.title = element_text(hjust = 0.5,size=15),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),axis.text.y= element_text(size=15,colour = 'black'),
        axis.text.x = element_text(size=15,colour = 'black'),legend.title = element_text('Antimicrobial peptide gene'),
        legend.text =element_text(size=11,colour = 'black') )
tight_plot
ggsave('tight_plot.png',tight_plot,width = 8,height = 8)


