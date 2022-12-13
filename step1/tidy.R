rm(list = ls()) 
options(stringsAsFactors = F) 
group1 <- read.table('../scz/PsychEncode_rnaseq_sczdata/together/output.txt')
clinical1 <- read.csv('../scz/PsychEncode_rnaseq_sczdata/PECCapstoneCollection_ClinicalData.csv')
clinical2 <- read.csv('../scz/PsychEncode_rnaseq_sczdata/PECCapstoneCollection_Tablesofsamples_RNAseq.csv')
group1_col <- colnames(group1)
library(stringr)
str_count(group1_col)
group1_coll <- str_sub(group1_col,start = 2,end = -26)
group1_coll
a <- str_replace(group1_coll,'\\.','-')
colnames(group1) <- a
match(a,clinical1$individualID)
clinical1$diagnosis[match(a,clinical1$individualID)]
group1_clinical <-clinical1$diagnosis[match(a,clinical1$individualID)] 
save(group1,group1_clinical,file = 'braingvex.Rdata')  
######################
group2 <- read.table('../scz/PsychEncode_rnaseq_sczdata/together/CMC_genecounts.txt')
group2_col <- colnames(group2)
match(group2_col,clinical2$name)
clinical2$individualID[match(group2_col,clinical2$name)]
match(clinical2$individualID[match(group2_col,clinical2$name)],clinical1$individualID)
clinical1$diagnosis[match(clinical2$individualID[match(group2_col,clinical2$name)],clinical1$individualID)]
group2_clinical <- clinical1$diagnosis[match(clinical2$individualID[match(group2_col,clinical2$name)],clinical1$individualID)]
save(group2,group2_clinical,file = 'CMC_genecount.Rdata')
#############################################################
group3 <- read.table('../scz/PsychEncode_rnaseq_sczdata/together/HBCC.txt')
group3_col <- colnames(group3)
match(group3_col,clinical2$name)
clinical2$individualID[match(group3_col,clinical2$name)]
match(clinical2$individualID[match(group3_col,clinical2$name)],clinical1$individualID)
clinical1$diagnosis[match(clinical2$individualID[match(group3_col,clinical2$name)],clinical1$individualID)]
group3_clinical <- clinical1$diagnosis[match(clinical2$individualID[match(group3_col,clinical2$name)],clinical1$individualID)]
save(group3,group3_clinical,file = 'HBCC_genecount.Rdata')
################################################################
group4 <- read.table('../scz/PsychEncode_rnaseq_sczdata/together/LIBD_out.txt',header = TRUE)
group4_col <- colnames(group4)
match(group4_col,clinical2$name)
clinical2$individualID[match(group4_col,clinical2$name)]
match(clinical2$individualID[match(group4_col,clinical2$name)],clinical1$individualID)
clinical1$diagnosis[match(clinical2$individualID[match(group4_col,clinical2$name)],clinical1$individualID)]
group4_clinical <- clinical1$diagnosis[match(clinical2$individualID[match(group4_col,clinical2$name)],clinical1$individualID)]
save(group4,group4_clinical,file = 'LIBD_genecount.Rdata')
