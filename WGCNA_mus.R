install.packages('WGCNA')#安装WGCNA包
library(WGCNA)#加载WGCNA包
library(stringr)
#打开多线程 
enableWGCNAThreads(6)
expression_all<-read.delim('Expression.txt',row.names = 1, sep = '\t',stringsAsFactors = F, 
                           check.names = F)
datExpr<-expression_all[,c(1:3,7:9,16:18,25:27,28:30)]
datExpr_1<-datExpr[which(rowSums(datExpr) > 1),]
datExpr_1<-t(datExpr_1[order(apply((datExpr_1),1,mad),decreasing = T),])
samplename<-as.data.frame(row.names(datExpr_1))
datTraits<-read.delim('detail.txt',sep = '\t',stringsAsFactors = T,check.names = F)
row.names(datTraits)<-datTraits$Samplename
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
soft <- pickSoftThreshold(datExpr_1, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(soft$fitIndices[,1], -sign(soft$fitIndices[,3])*soft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(soft$fitIndices[,1], soft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(soft$fitIndices[,1], soft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##corr matrix Modules
nGenes = ncol(datExpr_1)
nSamples = nrow(datExpr_1)
net<-blockwiseModules(datExpr_1,power = soft$powerEstimate,maxBlockSize = 10000,
                      TOMType = 'unsigned',minModuleSize = 30,reassignThresholdPS = 0,
                      mergeCutHeight = 0.25,numericLabels = TRUE,pamRespectsDendro = FALSE,
                      savelTOMs = T,verbose = 3,corType = 'pearson',loadTOMs=T)
table(net$colors)
#plotting the color label
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrograms and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, 
                    guideHang = 0.05)
#net heatmap
nGenes = ncol(datExpr_1)
nSamples = nrow(datExpr_1)
geneTree = net$dendrograms[[1]]
dissTOM <- 1-(as.matrix(TOMsimilarityFromExpr(datExpr_1, power = 14)))#一条命令到dissTOM
#TOM <- TOMsimilarityFromExpr(datExpr_1, power = 14)
#TOM<-as.matrix(TOM)
#dissTOM<- 1-TOM
#plotTOM = dissTOM^7
diag(plotTOM) = NA
moduleColors <- labels2colors(net$colors)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")#所有基因的TOM图
nSelect = 6000
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
sizeGrWindow(10,10)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
library(gplots)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
#sample hculust
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(datExpr_1), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
sample_colors <- numbers2colors(as.numeric(factor(datTraits$Days)), 
                                colors = c("white","blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

# Quantifying module–trait associations
# calculate module Eigengenes in all samples(ME矩阵)
design<-model.matrix(~0+ datTraits$Days)
colnames(design)=levels(datTraits$Days)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_1, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed (50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

design_1<-model.matrix(~0+ datTraits$Treat)
colnames(design_1)=levels(datTraits$Treat)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0_1 = moduleEigengenes(datExpr_1, moduleColors)$eigengenes
MEs_1 = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor_1 = cor(MEs, design_1 , use = "p");
moduleTraitPvalue_1 = corPvalueStudent(moduleTraitCor_1, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix_1 = paste(signif(moduleTraitCor_1, 2), "\n(",
                   signif(moduleTraitPvalue_1, 1), ")", sep = "");
dim(textMatrix_1) = dim(moduleTraitCor_1)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor_1,
               xLabels = colnames(design_1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = blueWhiteRed (50),
               textMatrix = textMatrix_1,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#days30, darkgrey color
#names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr_1, MEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
## 只有连续型性状才能只有计算
## 这里把是否属于 days30 表型这个变量用0,1进行数值化。
days30 = as.data.frame(design[,4])
names(days30) = "days30"
geneTraitSignificance = as.data.frame(cor(datExpr_1, days30, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(days30), sep="")
names(GSPvalue) = paste("p.GS.", names(days30), sep="")
#把两个相关性矩阵联合起来,指定感兴趣模块进行分析,选取饲喂时间正相关性
module = "royalblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for days30",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,abline = T,
                   lmFnc = lm)
#days30, tan color
#names (colors) of the modules
#modNames = substring(names(MEs), 3)
#geneModuleMembership = as.data.frame(cor(datExpr_1, MEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
#MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#names(geneModuleMembership) = paste("MM", modNames, sep="")
#names(MMPvalue) = paste("p.MM", modNames, sep="")
## 只有连续型性状才能只有计算
## 这里把是否属于 days20 表型这个变量用0,1进行数值化。
#days30 = as.data.frame(design[,4])
#names(days30) = "days30"
#geneTraitSignificance = as.data.frame(cor(datExpr_1, days30, use = "p"))
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
#names(geneTraitSignificance) = paste("GS.", names(days30), sep="")
#names(GSPvalue) = paste("p.GS.", names(days30), sep="")
#把两个相关性矩阵联合起来,指定感兴趣模块进行分析
#module = "tan"
#column = match(module, modNames);
#moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                   abs(geneTraitSignificance[moduleGenes, 1]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = "Gene significance for days30",
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,abline = T,
#                   lmFnc = lm)

#模块中基因热图与表达bar—plot
sizeGrWindow(8,7)
which.module = 'darkgrey'
ME = MEs[,paste('ME',which.module,sep = '')]
par(mfrow = c(2,1),mar=c(0.3,5.5,3,2))
plotMat(t(scale(datExpr_1[,moduleColors == which.module])),nrgcols = 30,rlabels = F,
        rcols = which.module,main = which.module,cex.main = 2)
par(mar=c(5,4.2,0,0.7))
barplot(ME,col = which.module, main = '',cex.main = 2,ylab = 'eigengene expression',xlab = 'sample')

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr_1, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属于 days30 表型这个变量用0,1进行数值化。
days30 = as.data.frame(design[,4])
names(days30) = "days30"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, days30))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle= 90)
# Plot the dendrogram
sizeGrWindow(6,6)
par(cex = 1.0)
## 模块的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

#提取royalblue模块基因
module = "darkgrey"
# Select module probes
probes = colnames(datExpr_1) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
head(modProbes)#之后可以进行GO与KEGG分析
#提取tan模块基因
module = "turquoise"
#Select module probes
probes = colnames(datExpr_1) 
inModule = (moduleColors==module)
modProbes_1 = probes[inModule]
head(modProbes_1)#之后可以进行GO与KEGG分析
#probe gene GO and KEGG analysis_royalblue
library(clusterProfiler)
library(org.Mm.eg.db)
GO_gene<- enrichGO(gene = modProbes,OrgDb= org.Mm.eg.db,
                         keyType= 'ENSEMBL',ont= "ALL",pAdjustMethod = "BH",pvalueCutoff  = 1)
GO_gene_detail<-as.data.frame(GO_gene)
GO_gene_BP_detail<-GO_gene_detail[which(GO_gene_detail$ONTOLOGY == 'BP'),]
GO_gene_CC_detail<-GO_gene_detail[which(GO_gene_detail$ONTOLOGY == 'CC'),]
GO_gene_MF_detail<-GO_gene_detail[which(GO_gene_detail$ONTOLOGY == 'MF'),]
GO_gene_BP_detail['total']<-45
GO_gene_MF_detail['total']<-44
GO_gene_CC_detail['total']<-46
GO_gene_BP_detail['Richfactor']<-(GO_gene_BP_detail$Count)/(GO_gene_BP_detail$total)
GO_gene_CC_detail['Richfactor']<-(GO_gene_CC_detail$Count)/(GO_gene_CC_detail$total)
GO_gene_MF_detail['Richfactor']<-(GO_gene_MF_detail$Count)/(GO_gene_MF_detail$total)
GO_gene_BP_detail<-GO_gene_BP_detail[order(GO_gene_BP_detail$Count,decreasing = T),]
GO_gene_MF_detail<-GO_gene_MF_detail[order(GO_gene_MF_detail$Count,decreasing = T),]
GO_gene_CC_detail<-GO_gene_CC_detail[order(GO_gene_CC_detail$Count,decreasing = T),]
library(ggplot2)
GO_royalblue_gene_BP_detail_top20plot<-ggplot(GO_royalblue_gene_BP_detail[c(1:20),],aes(Richfactor,Description))+
        geom_point(aes(size=Count,color= pvalue))+scale_color_gradient(low = "red",high = "blue")+
        labs(title="Biological Process Top20")+theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
              axis.text.x = element_text(size=12,colour = 'black'))
GO_royalblue_gene_BP_detail_top20plot
ggsave('GO_royalblue_gene_BP_detail_top20plot.tiff',GO_royalblue_gene_BP_detail_top20plot,width = 10,height = 8)
GO_royalblue_gene_MF_detail_top20plot<-ggplot(GO_royalblue_gene_MF_detail[c(1:20),],aes(Richfactor,Description))+
        geom_point(aes(size=Count,color= pvalue))+scale_color_gradient(low = "red",high = "blue")+
        labs(title="Molecular Function Top20")+theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
              axis.text.x = element_text(size=12,colour = 'black'))
GO_royalblue_gene_MF_detail_top20plot
ggsave('GO_royalblue_gene_MF_detail_top20plot.tiff',GO_royalblue_gene_MF_detail_top20plot,width = 8,height = 8)
GO_royalblue_gene_CC_detail_top20plot<-ggplot(GO_royalblue_gene_CC_detail[c(1:20),],aes(Richfactor,Description))+
        geom_point(aes(size=Count,color= pvalue))+scale_color_gradient(low = "red",high = "blue")+
        labs(title="Cellular Component Top20")+theme(plot.title = element_text(hjust = 0.5,size=15),
                                                     axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
                                                     axis.text.x = element_text(size=12,colour = 'black'))
GO_royalblue_gene_CC_detail_top20plot
ggsave('GO_royalblue_gene_CC_detail_top20plot.tiff',GO_royalblue_gene_CC_detail_top20plot,width = 8,height = 8)
KEGG_green_geneid<-bitr(unique(modProbes), fromType = "ENSEMBL",toType = c( "ENTREZID"),
                        OrgDb = org.Mm.eg.db)
KEGG_green_gene_KEGG<- enrichKEGG(gene = KEGG_green_geneid$ENTREZID,organism = "mouse",
                              keyType = "kegg",pvalueCutoff = 1,qvalueCutoff = 1,
                              pAdjustMethod = "BH")
KEGG_green_gene_detail<-as.data.frame(KEGG_green_gene_KEGG)
KEGG_green_gene_detail['total']<-435
KEGG_green_gene_detail['Richfactor']<-(KEGG_green_gene_detail$Count)/(KEGG_green_gene_detail$total)
KEGG_green_gene_detail<-KEGG_green_gene_detail[order(KEGG_green_gene_detail$Count,decreasing = T),]
KEGG_green_gene_detail_top20plot<-ggplot(KEGG_green_gene_detail[c(1:20),],aes(Richfactor,Description))+
        geom_point(aes(size=Count,color= pvalue))+scale_color_gradient(low = "red",high = "blue")+labs(title="KEGG Pathway Top20")+
        theme(plot.title = element_text(hjust = 0.5,size=15),
              axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),axis.text.y= element_text(size=12,colour = 'black'),
              axis.text.x = element_text(size=12,colour = 'black'))
KEGG_green_gene_detail_top20plot

pathway<-grep('pathway',KEGG_green_gene_detail$Description)
KEGG_green_gene_pathway<-KEGG_green_gene_detail[pathway,]

#模块导出
# Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExpr_1, power = soft$powerEstimate)
# Select module
module <- "darkgrey"
# Select module probes
probes = colnames(datExpr_1) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#export to cytoscape
cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule])
module <- "turquoise"
# Select module probes
probes = colnames(datExpr_1) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#export to cytoscape
cyt_1 = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule])

#Hub gene
connectivity<-abs(cor(datExpr_1,use="p"))^soft #跟TOM一样
Alldegrees=intramodularConnectivity(connectivity, mergedColors) #计算每个基因module内的连接度和module外的连接度以及total连接度。
Alldegrees$gene=rownames(Alldegrees)
datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")#计算MM。
GS=cor(datExpr,design,use="p") #注意，如果是design是离散变量，一定要将design放在后面
combin_data=cbind(Alldegrees,datKME,GS1)
write.csv(combin_data,"combine_hub.csv")

#
Top_hubgene<-chooseTopHubInEachModule(datExpr_1,mergedColors,omitColors = "grey",power = 14,type = "signed")


#
library(clusterProfiler)
library(org.Mm.eg.db)
cytoscape_gene<-bitr(unique(modProbes), fromType = "ENSEMBL",toType = c( "SYMBOL"),
                              OrgDb = org.Mm.eg.db)
cytoscape_gene_edge<-read.delim('CytoscapeInput-edges-darkgrey.txt',sep = '\t',check.names = F,
                                stringsAsFactors = F)
cytoscape_gene['fromNode']<-cytoscape_gene$ENSEMBL
cytoscape_gene_edge<-merge(cytoscape_gene_edge,cytoscape_gene,by='fromNode',all.x=T)
cytoscape_gene_edge['From']<-cytoscape_gene_edge$SYMBOL
cytoscape_gene_edge<-cytoscape_gene_edge[,-c(1,7:8)]
cytoscape_gene['toNode']<-cytoscape_gene$ENSEMBL
cytoscape_gene_edge<-merge(cytoscape_gene_edge,cytoscape_gene,by='toNode',all.x=T)
cytoscape_gene_edge['To']<-cytoscape_gene_edge$SYMBOL
cytoscape_gene_edge<-cytoscape_gene_edge[,-c(1,7:9)]
write.csv(file = 'cytoscape_gene_edge.csv',cytoscape_gene_edge)

##
D10NTC_D30TOA_DEG_Syembol<-read.csv('D10NTC_D30TOA_DEG_1.csv')
D10NTC_D30TOA_DEG_Syembol<-D10NTC_D30TOA_DEG_Syembol[,c(8:9)]
cytoscape_gene_edge['SYMBOL']<-cytoscape_gene_edge$From
cytoscape_gene_edge<-cytoscape_gene_edge[,-c(3:4)]
cytoscape_gene_edge<-merge(cytoscape_gene_edge,D10NTC_D30TOA_DEG_Syembol,by = 'SYMBOL',all.x =T )
cytoscape_gene_edge<-cytoscape_gene_edge[,-1]
write.csv(file = 'cytoscape_gene_edge.csv',cytoscape_gene_edge)
node<-read.csv('cytoscape_gene_edge.csv default node.csv')
D10NTC_D30TOA_DEG_Syembol['name']<-D10NTC_D30TOA_DEG_Syembol$SYMBOL
node<-merge(node,D10NTC_D30TOA_DEG_Syembol,by = 'name',all.x = T)
write.csv(file = 'cytoscape_gene_default node.csv',node)
##