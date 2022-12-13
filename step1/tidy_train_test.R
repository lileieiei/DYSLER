rm(list=ls())
load("HBCC_genecount.Rdata") # group3_ TEST_SET !
load('LIBD_genecount.Rdata') # group4
load('CMC_genecount.Rdata') # group2
load('braingvex.Rdata')  #group1
#table(group3_clinical)
#table(group1_clinical)
#identical(rownames(group1),rownames(group2))
#identical(rownames(group2),rownames(group4))
ss <- cbind(group1,group2)
ss <- cbind(ss,group4)
ss_clin <- c(group1_clinical,group2_clinical,group4_clinical)
table(ss_clin)
aa <- unique(ss_clin)
ss_clin
ss_clin[which(ss_clin%in%aa[c(2,4,5,6)])]
clin_keep <- ss_clin[-which(ss_clin%in%aa[c(2,4,5,6)])] #keep data
group_keep <- ss[,-which(ss_clin%in%aa[c(2,4,5,6)])]  #keep data
table(clin_keep)
table(ss_clin)
table(group3_clinical)
which(group3_clinical%in%aa[c(2,4,5,6)])
HBCC_exp <- group3[,-which(group3_clinical%in%aa[c(2,4,5,6)])]
HBCC_clin <- group3_clinical[-which(group3_clinical%in%aa[c(2,4,5,6)])]
save(HBCC_exp,HBCC_clin,group_keep,clin_keep,file="final_keep.Rdata")
