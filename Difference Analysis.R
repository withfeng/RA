## 导入R包
library(limma)
library(dplyr)

df <- read.table("Counts_data.txt", header = T, sep = "\t", row.names = 1, check.names = F)
head(df)

list <- c(rep("Treat", 66), rep("CK",148)) %>% factor(., levels = c("CK", "Treat"), ordered = F)
##--------------
> head(list)
  CK Treat
1  0     1
2  0     1
3  0     1
4  0     1
5  0     1
6  0     1
list <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
colnames(list) <- c("CK", "Treat")
df.fit <- lmFit(df, list)  ## 数据与list进行匹配

df.matrix <- makeContrasts(tumor - normal, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)

## 导出所有的差异结果
nrDEG = na.omit(tempOutput) ## 去掉数据中有NA的行或列
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut.csv")  
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")

## 我们使用|logFC| > 1，padj < 0.05（矫正后P值）
foldChange = 1
padj = 0.05
## 筛选出所有差异基因的结果
All_diffSig <- diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
#---------------------
> dim(All_diffSig)

diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup.csv")
#
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig < -foldChange)),]
write.csv(diffdown, "diffdown.csv")

## 导入R包
library(ggplot2)
library(ggrepel)
##  绘制火山图
## 进行分类别
logFC <- diffsig$logFC
deg.padj <- diffsig$P.Value
data <- data.frame(logFC = logFC, padj = deg.padj)
data$group[(data$padj > 0.05 | data$padj == "NA") | (data$logFC < foldChange) & data$logFC > -foldChange] <- "Not"
data$group[(data$padj <= 0.05 & data$logFC > 1)] <-  "Up"
data$group[(data$padj <= 0.05 & data$logFC < -1)] <- "Down"
x_lim <- max(logFC,-logFC)

# 开始绘图
pdf('volcano.pdf',width = 7,height = 6.5)  ## 输出文件
label = subset(diffsig,P.Value <0.05 & abs(logFC) > 0.5)
label1 = rownames(label)

colnames(diffsig)[1] = 'log2FC'
Significant=ifelse((diffsig$P.Value < 0.05 & abs(diffsig$log2FC)> 0.5), ifelse(diffsig$log2FC > 0.5,"Up","Down"), "Not")

ggplot(diffsig, aes(log2FC, -log10(P.Value)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#0072B5","grey","#BC3C28"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-0.5,0.5), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")+
  theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
  str(diffsig, max.level = c(-1, 1))+theme_bw()

dev.off()

## 导入R包
library(pheatmap)

## 提取差异基因的表达量
DEG_id <- read.csv("all.diffsig.csv", header = T)  # 读取差异基因的文件
head(DEG_id)
## 匹配差异基因的表达量
DEG_id <- unique(DEG_id$X)
DEG_exp <- df[DEG_id,]
hmexp <- na.omit(DEG_exp)

## 样本注释信息 
annotation_col <- data.frame(Group = factor(c(rep("Treat", 66), rep("CK",148))))
rownames(annotation_col) <- colnames(hmexp)

##  绘制热图 
pdf(file = "heatmap02.pdf", height = 8, width = 12)
pheatmap(hmexp,
              annotation_col = annotation_col,
              color = colorRampPalette(c("green","black","red"))(50),
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
print(p)
dev.off()
