rm(list=ls())
library(tinyarray)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(stringi)##加载包
cancer_type="SCZ"
load(paste0(cancer_type,"DEG.Rdata"))
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
#table(!rownames(dd3_up)%in%rownames(dd1_up))
dd2_up_add <- dd2_up[!rownames(dd2_up)%in%rownames(dd1_up),]
dd1_up <- rbind(dd1_up,dd2_up_add)
dd3_up_add <- dd3_up[!rownames(dd3_up)%in%rownames(dd1_up),]
dd1_up <- rbind(dd1_up,dd3_up_add)
dd2_down_add <- dd2_down[!rownames(dd2_down)%in%rownames(dd1_down),]
dd1_down <- rbind(dd1_down,dd2_down_add)
dd3_down_add <- dd3_down[!rownames(dd3_down)%in%rownames(dd1_down),]
dd1_down <- rbind(dd1_down,dd3_down_add)
gene_up <- dd1_up  #2312 个
gene_down <- dd1_down #916个
gene_up <- stri_sub(rownames(dd1_up),1,15)
gene_down <- stri_sub(rownames(dd1_down),1,15)
gene_up_entrezid <- bitr(gene_up,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)[,2]
gene_down_entrezid <- bitr(gene_down,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)[,2]
save(gene_up,gene_down,file = "degfor_keggandgo.Rdata")
if(F){
  #CC
  ego_CC_up <- enrichGO(gene = gene_up_entrezid,
                      OrgDb= org.Hs.eg.db,
                      ont = "CC",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  #BP
  ego_BP_up <- enrichGO(gene = gene_up_entrezid,
                      OrgDb= org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  #MF
  ego_MF_up <- enrichGO(gene = gene_up_entrezid,
                      OrgDb= org.Hs.eg.db,
                      ont = "MF",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = TRUE)
  save(gene_up_entrezid,ego_CC_up,ego_BP_up,ego_MF_up,file = "ego_SCZ_up.Rdata")
}

#低表达基因
if(F){
  #CC
  ego_CC_down <- enrichGO(gene = gene_down_entrezid,
                        OrgDb= org.Hs.eg.db,
                        ont = "CC",
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  #BP
  ego_BP_down <- enrichGO(gene = gene_down_entrezid,
                        OrgDb= org.Hs.eg.db,
                        ont = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  #MF
  ego_MF_down <- enrichGO(gene = gene_down_entrezid,
                        OrgDb= org.Hs.eg.db,
                        ont = "MF",
                        pAdjustMethod = "BH",
                        minGSSize = 1,
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  save(gene_down_entrezid,ego_CC_down,ego_BP_down,ego_MF_down,file = "ego_SCZ_down.Rdata")
}
#高表达基因kegg
kk.diff_up <- enrichKEGG(gene         = gene_up_entrezid,
                       organism     = 'hsa')
#低表达基因kegg
kk.diff_down <- enrichKEGG(gene         = gene_down_entrezid,
                         organism     = 'hsa')
kegg_up <- as.data.frame(kk.diff_up)
kegg_down <- as.data.frame(kk.diff_down)
bp_up <- as.data.frame(ego_BP_up)
bp_down <- as.data.frame(ego_BP_down)
library(stringr)
bp_up_target <- unique(unlist(str_split(bp_up[,'geneID'],'/')))
bp_down_target <- unique(unlist(str_split(bp_down[,'geneID'],'/')))
kegg_up_target <- unique(unlist(str_split(kegg_up[,'geneID'],'/')))
kegg_down_target <- unique(unlist(str_split(kegg_down[,'geneID'],'/')))
kegg_up_target=bitr(kegg_up_target,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)[,2]
kegg_down_target=bitr(kegg_down_target,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db)[,2]
up_all=unique(c(bp_up_target,kegg_up_target))
down_all=unique(c(bp_down_target,kegg_down_target))
############################
load("ego_SCZ_up.Rdata")
load("ego_SCZ_down.Rdata")


###导出文件
write.csv(kegg_up,file='up_regulated_kegg.csv')
write.csv(kegg_down,file='down_regulated_kegg.csv')
write.csv(bp_down,file='down_regulated_bp.csv')
write.csv(bp_up,file='up_regulated_bp.csv')
