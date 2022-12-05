#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#���ð�
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        #pֵ��������
qvalueFilter=0.05        #�������pֵ��������

setwd("C:\\Users\\lexb4\\Desktop\\geoSeq\\14.KEGG")        #���ù���Ŀ¼
rt=read.table("id.txt",sep="\t",header=T,check.names=F)    #��ȡid.txt�ļ�
rt=rt[is.na(rt[,"entrezID"])==F,]                          #ȥ������idΪNA�Ļ���
colnames(rt)[1]="Gene"
gene=rt$entrezID
geneFC=rt$logFC
#geneFC=2^geneFC
names(geneFC)=gene

#������ɫ����
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#kegg��������
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$Gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#���渻�����
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#������ʾTerm��Ŀ
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#��״ͼ
pdf(file="barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)
dev.off()

#����ͼ
pdf(file="bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
dev.off()

#Ȧͼ
pdf(file="circos.pdf",width = 11,height = 7)
kkx=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kkx, foldChange=geneFC,showCategory = 5, circular = TRUE, colorEdge = TRUE,node_label="all")
dev.off()


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio