rm(list=ls())
kegg_up <- read.csv('up_regulated_kegg.csv')
kegg_down <- read.csv('down_regulated_kegg.csv',row.names = 1)
bp_down <- read.csv('down_regulated_bp.csv',row.names = 1)
bp_up <- read.csv('up_regulated_bp.csv',row.names = 1)
kegg_up_num <-c(2,4,5,6,10,11,13)-1 
kegg_down_num <- 1
bp_up_num <- c(3,4,5,9,11,17,28,30,32,33,41,45,47,51,53,55,56,58,59,60)-1
bp_down_num <- c(16,24,37,41,44,47,53,55,74,75,82)-1
kegg_up_shift <- kegg_up[kegg_up_num,]
kegg_down_shift <- kegg_down[kegg_down_num,]
bp_up_shift <- bp_up[bp_up_num,]
bp_down_shift <- bp_down[bp_down_num,]
#####
library(stringr)
bp_up_target <- unique(unlist(str_split(bp_up_shift[,'geneID'],'/')))
bp_down_target <- unique(unlist(str_split(bp_down_shift[,'geneID'],'/')))
kegg_up_target <- unique(unlist(str_split(kegg_up_shift[,'geneID'],'/')))
kegg_down_target <- unique(unlist(str_split(kegg_down_shift[,'geneID'],'/')))
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
kegg_up_target <- bitr(kegg_up_target,fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb =org.Hs.eg.db)[,2]
kegg_down_target <- bitr(kegg_down_target,fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb =org.Hs.eg.db)[,2]
library(tinyarray)

x=list(kegg_up=kegg_up_target,kegg_down=kegg_down_target,bp_up=bp_up_target,bp_down=bp_down_target)
draw_venn(x,'tt')
symbol_all <- unique(c(kegg_up_target,kegg_down_target,bp_up_target,bp_down_target))
symbol=bitr(symbol_all,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb =org.Hs.eg.db)[,2]
write.table(symbol,file='tonglushaixuanjieguo.txt',row.names = FALSE,quote = FALSE)

rm(list=ls())
load('final_keep.Rdata')
library(stringr)
library(stringi)
ss <- stri_sub(rownames(group_keep),1,15)
ss_test <- stri_sub(rownames(HBCC_exp),1,15)
aa <- read.table('./tonglushaixuanjieguo.txt',header = T)[,1]
exp <- group_keep[na.omit(match(aa,ss)),]
exp_test <- HBCC_exp[na.omit(match(aa,ss_test)),]
expset <- as.data.frame(t(exp))
expset_test <- as.data.frame(t(exp_test))

expset$group=as.numeric(as.factor(clin_keep))-1
expset_test$group=as.numeric(as.factor(HBCC_clin))-1

sci_gene <- read.csv('./shift_fdr_science.csv',header = T)
ss_sci <- stri_sub(colnames(expset),1,15)[1:138]
shiftt <- ss_sci[ss_sci%in%sci_gene$ensembl_gene_id]
expsett=expset[,match(shiftt,ss_sci)]
expset_testt=expset_test[,match(shiftt,ss_sci)]
expsett$group=as.numeric(as.factor(clin_keep))-1
expset_testt$group=as.numeric(as.factor(HBCC_clin))-1

write.csv(expset_testt,'6_1data_exp_test_shifted.csv')
write.csv(expsett,'6_1data_exp_shifted.csv')
