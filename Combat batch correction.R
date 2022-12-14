library(limma)
library(sva)
library(tinyarray)

#文件名称
files=c("GSE141512.txt", "GSE59867.txt")      

geneList=list()
for(i in 1:length(files)){
	fileName=files[i]
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    header=unlist(strsplit(fileName, "\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
	fileName=files[i]
    header=unlist(strsplit(fileName, "\\.|\\-"))
    rt=read.table(fileName, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
    rt=avereps(data)
   colnames(rt)=paste0(header[1], "_", colnames(rt))
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
    	rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    if(i==1){
    	allTab=rt[intersectGenes,]
    }else{
    	allTab=cbind(allTab, rt[intersectGenes,])
    }
    
}

#整合数据
sample=read.table("sample.txt", header=F, sep="\t", check.names=F,row.names = 1)
allTab=allTab[,rownames(sample)]
allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file="raw.exp.txt", sep="\t", quote=F, col.names=F)

#batchType系列来源，modeType样本分组类型
batchType=as.character(as.data.frame(strsplit(rownames(sample), "_"))[1,])
modeType=as.vector(sample[,1])
mod=model.matrix(~as.factor(modeType))
normalizeTab=ComBat(allTab, batchType,mod, par.prior=TRUE)
normalizeTabout=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTabout, file="combat.exp.txt", sep="\t", quote=F, col.names=F)

#############矫正前
#样本分组
type=factor(x = colnames(allTab),labels =modeType)
summary(type)
pdf(file = "raw_PCA_Type.pdf",width=5.5,height = 4.5)
draw_pca(allTab,type)
dev.off()
#样本系列来源
type=factor(x = colnames(allTab),labels =batchType)
summary(type)
pdf(file = "raw_PCA_Group.pdf",width=5.5,height = 4.5)
draw_pca(allTab,type)
dev.off()
#############矫正后
#样本分组
type=factor(x = colnames(normalizeTab),labels =modeType)
summary(type)
pdf(file = "nor_PCA_Type.pdf",width=5.5,height = 4.5)
draw_pca(normalizeTab,type)
dev.off()
#样本系列来源
type=factor(x = colnames(normalizeTab),labels =batchType)
summary(type)
pdf(file = "nor_PCA_Group.pdf",width=5.5,height = 4.5)
draw_pca(normalizeTab,type)
dev.off()