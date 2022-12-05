######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("reshape2")


#���ð�
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

expFile="diffGeneExp.txt"           #��������
immFile="CIBERSORT-Results.txt"     #����ϸ���������ļ�
setwd("C:\\Users\\LENOVO\\Desktop\\cor-CIBRSORT")     #���ù���Ŀ¼

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#ȥ����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)

#��ȡ����ϸ������ļ����������ݽ�������
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]

#����Է���
outTab=data.frame()
for(cell in colnames(immune)){
	if(sd(immune[,cell])==0){next}
	for(gene in colnames(data)){
		x=as.numeric(immune[,cell])
		y=as.numeric(data[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
	}
}

#�����������ͼ
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=5)
ggplot(outTab, aes(Immune, Gene)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #ȥ������
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),   #x������
	      axis.text.y = element_text(size = 8, face = "bold")) +       #y������
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #����ͼ��
	scale_x_discrete(position = "bottom")      #X��������ʾλ��
dev.off()


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056