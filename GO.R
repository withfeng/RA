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

pvalueFilter=0.05         #pֵ��������
qvalueFilter=0.05         #�������pֵ��������

setwd("C:\\Users\\lexb4\\Desktop\\geoSeq\\13.GO")         #���ù���Ŀ¼
rt=read.table("id.txt",sep="\t",header=T,check.names=F)      #��ȡid.txt�ļ�
rt=rt[is.na(rt[,"entrezID"])==F,]                            #ȥ������idΪNA�Ļ���
gene=rt$entrezID
geneFC=rt$logFC
#geneFC=2^geneFC
names(geneFC)=gene

#������ɫ����
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#GO��������
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#���渻�����
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#������ʾTerm��Ŀ
showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

#��״ͼ
pdf(file="barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#����ͼ
pdf(file="bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
		
#Ȧͼ
pdf(file="circos.pdf",width = 9,height = 6.5)
cnet=cnetplot(kk, foldChange=geneFC, showCategory = 5, circular = TRUE, colorEdge = TRUE)
print(cnet)
dev.off()


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio