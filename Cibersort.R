library(psych)
library(ggcorrplot)

setwd("C:\\Users\\LENOVO\\Desktop\\666")

# data1=read.table("result.xls",header = T,sep = "\t",check.names =F)
# data2=read.table("123.xls",header = T,sep = "\t",check.names =F)
data1 = readxl::read_xlsx("C:\\Users\\LENOVO\\Desktop\\666\\result.xls",sheet = 1)
data2 = readxl::read_xlsx("C:\\Users\\LENOVO\\Desktop\\666\\123.xls",sheet = 1)

data1 = as.data.frame(data1)
rownames(data1) = data1$id
data1 = data1[,-1]
data2 = as.data.frame(data2)
rownames(data2) = data2$id
data2 = data2[,-1]

inter_sample = intersect(rownames(data1),colnames(data2))
data1 = data1[inter_sample,]
data2 = data2[,inter_sample]

# rownames(data2)=data2$id
mygene=c("CCR5","CD27","IL6", "LCK", "SELL","CD2","GZMB","IL2RG","IL7R","CD40LG","CXCL10","IL15","ZAP70")
data3=t(data2[mygene,])
# rownames(data3)=NULL
# data3=as.matrix(data3)
# data3=t(data3)
# colnames(data3)=data3[1,]
# data3=data3[-1,]
# class(data3)
# data1=t(data1)
# colnames(data1)=data1[1,]
# data1=data1[-1,]
# data1=apply(data1,2,as.numeric)
# data3=apply(data3,2,as.numeric)

mycor<- corr.test(data3,data1,method = 'spearman')
rvalue <- mycor$r
pvalue <- mycor$p

p <- ggcorrplot(t(rvalue), show.legend = T, 
                p.mat = t(mycor$p.adj), digits = 2,  sig.level = 0.05,insig = 'blank',lab = T)
ggsave("1.pdf", p, width = 20, height = 8)

