#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        #p值过滤条件
qvalueFilter=0.05        #矫正后的p值过滤条件

setwd("C:\\Users\\lexb4\\Desktop\\geoSeq\\14.KEGG")        #设置工作目录
rt=read.table("id.txt",sep="\t",header=T,check.names=F)    #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                          #去除基因id为NA的基因
colnames(rt)[1]="Gene"
gene=rt$entrezID
geneFC=rt$logFC
#geneFC=2^geneFC
names(geneFC)=gene

#定义颜色类型
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#定义显示Term数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

#圈图
pdf(file="circos.pdf",width = 11,height = 7)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, foldChange=geneFC,showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
