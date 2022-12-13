rm(list=ls())
load('final_keep.Rdata')
#table(clin_keep)
#table(HBCC_clin)
#library(devtools)
#devtools::install_github("xjsun1221/tinyarray")
library(tinyarray)
#table(complete.cases(group_keep))
exp <- round(group_keep) # data

new_exp = exp[apply(exp, 1, function(x) sum(x > 1) > 50), ]


group_list = factor(clin_keep,levels = c("Control","Schizophrenia"))
library(DESeq2)
cancer_type="SCZ"
colData <- data.frame(row.names =colnames(new_exp), 
                      condition=group_list)
if(!file.exists(paste0(cancer_type,"dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = new_exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(cancer_type,"dd.Rdata"))
}
load(paste0(cancer_type,"dd.Rdata"))

# compare
res <- results(dds, contrast = c("condition",rev(levels(group_list))))
resOrdered <- res[order(res$pvalue),] # 
DEG <- as.data.frame(resOrdered)
head(DEG)


logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
#logFC_cutoff <- 0.5
k1 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
k2 = (DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)

DESeq2_DEG <- DEG
a <- with(DESeq2_DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
#####
# ss <- t(new_exp["5221",])
# tr <- data.frame(value=ss,group=group_list)
# library(ggplot2)
# ggplot(data=tr,aes(x=group,y=X5221,fill=group)) + geom_boxplot()
########
library(edgeR)

dge <- DGEList(counts=new_exp,group=group_list)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 

design <- model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)

dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast=c(-1,1)) 

DEG=topTags(fit2, n=nrow(exp))
DEG=as.data.frame(DEG)
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

#logFC_cutoff <- 0.5
k1 = (DEG$PValue < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$PValue < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

head(DEG)
table(DEG$change)
edgeR_DEG <- DEG
b <- with(edgeR_DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
###limma----
library(limma)

design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(new_exp)

dge <- DGEList(counts=new_exp)
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(group_list)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)

logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )

#logFC_cutoff <- 0.5
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)

limma_voom_DEG <- DEG
c <- with(limma_voom_DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
tj = data.frame(deseq2 = as.integer(table(DESeq2_DEG$change)),
                edgeR = as.integer(table(edgeR_DEG$change)),
                limma_voom = as.integer(table(limma_voom_DEG$change)),
                row.names = c("down","not","up")
);tj
save(DESeq2_DEG,edgeR_DEG,limma_voom_DEG,group_list,tj,file = paste0(cancer_type,"DEG.Rdata"))
exp <- new_exp
library(tinyarray)
exp[1:4,1:4]
dat = log2(exp+1)
pca.plot = draw_pca(dat,group_list);pca.plot
save(pca.plot,file = paste0(cancer_type,"pcaplot.Rdata"))

cg1 = rownames(DESeq2_DEG)[DESeq2_DEG$change !="NOT"]
cg2 = rownames(edgeR_DEG)[edgeR_DEG$change !="NOT"]
cg3 = rownames(limma_voom_DEG)[limma_voom_DEG$change !="NOT"]

h1 = draw_heatmap(dat[cg1,],group_list,legend = T,annotation_legend=T)
h2 = draw_heatmap(dat[cg2,],group_list,legend = T,annotation_legend=T)
h3 = draw_heatmap(dat[cg3,],group_list,legend = T,annotation_legend=T)

v1 = draw_volcano(DESeq2_DEG,pkg = 1,logFC_cutoff =a)
v2 = draw_volcano(edgeR_DEG,pkg = 2,logFC_cutoff = b)
v3 = draw_volcano(limma_voom_DEG,pkg = 3,logFC_cutoff =c)

library(patchwork)
library(ggplot2)
(h1 + h2 + h3) / (v1 + v2 + v3)+plot_layout(guides = 'collect') 
#&theme(legend.position = "none")
ggsave(paste0(cancer_type,"heat_vo_91.png"),width = 20,height = 10,dpi=900)
###############################
#################### veen
load("SCZDEG.Rdata")
dd1 <- DESeq2_DEG[DESeq2_DEG$change!="NOT",]
dd2 <- edgeR_DEG[edgeR_DEG$change!="NOT",]
dd3 <- limma_voom_DEG[limma_voom_DEG$change!="NOT",]
dd1 <- dd1[,c(2,7)]
dd2 <- dd2[,c(1,6)]
dd3 <- dd3[,c(1,7)]
colnames(dd1) <- c("logFC","change")
dd1_up <- dd1[dd1$change=="UP",]
dd2_up <- dd2[dd2$change=="UP",]
dd3_up <- dd3[dd3$change=="UP",]
dd1_down <- dd1[dd1$change=="DOWN",]
dd2_down <- dd2[dd2$change=="DOWN",]
dd3_down <- dd3[dd3$change=="DOWN",]
library(tinyarray)
x=list(DESeq2_up=rownames(dd1_up),edgeR_up=rownames(dd2_up),limma_voom_up=rownames(dd3_up))

library(ggplot2)
draw_venn(x,'up_regulated_biomarker')
ggsave('up_regulated_biomarker.png',dpi=600)
dev.off()


x=list(DESeq2_down=rownames(dd1_down),edgeR_down=rownames(dd2_down),limma_voom_down=rownames(dd3_down))

draw_venn(x,'down_regulated_biomarker')
ggsave('down_regulated_biomarker.png',dpi=600)
dev.off()
