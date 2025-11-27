##colocalization
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(data.table)

#input full summary data of ENOC
data<-fread("chr6_ENOC.txt",header = T)

#step1 eQTL 
#rs880347 

gwas1 <- fread("eqtl_pancreas.txt",header = T)
gwas1 <- gwas1[gwas1$chr == "chr6",]
gwas1 <- gwas1[gwas1$pos> posx-500000   &  gwas1$pos < posx+500000,]
gwas1$samplesize<-round(gwas1$mac/gwas1$maf/2)
gwas1<-data.frame(gwas1$snp,
                  gwas1$pos,
                  gwas1$beta,
                  gwas1$se,
                  gwas1$eaf,
                  gwas1$p,
                  gwas1$ref,
                  gwas1$alt,
                  gwas1$phenotype_id,
                  gwas1$samplesize
)
colnames(gwas1)<-c("SNP","pos","beta","se","eaf","pval","other_allele","effect_allele","gene_id","N")
gwas2<-data
gwas2<- gwas2[gwas2$Position > posy-500000  &  gwas2$Position < posy+500000,]

colnames(gwas2)<-c("SNP","pos","beta","se","eaf","pval","other_allele","effect_allele","N")
dat_merge <- merge(gwas1,gwas2,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))

dat_merge$label_a1 <- ifelse(dat_merge$effect_allele_gwas1 == dat_merge$effect_allele_gwas2,"yes","no")
dat_merge$label_a1a2 <- ifelse(dat_merge$effect_allele_gwas1 == dat_merge$other_allele_gwas2,"yes","no")
dat_merge$label_a2a1 <- ifelse(dat_merge$other_allele_gwas1 == dat_merge$effect_allele_gwas2,"yes","no")
dat_merge$effect_allele_gwas2 <- ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                          dat_merge$effect_allele_gwas1,
                                          dat_merge$effect_allele_gwas2 )
dat_merge$other_allele_gwas2 <- ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                         dat_merge$other_allele_gwas1,
                                         dat_merge$other_allele_gwas2 )
dat_merge$beta_gwas2 <-ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                -dat_merge$beta_gwas2 ,
                                dat_merge$beta_gwas2)
dat_merge[,c("label_a1","label_a2","label_a1a2","label_a2a1")] <- NULL

dat_merge$MAF_gwas1   <- ifelse(dat_merge$eaf_gwas1   < 0.5,dat_merge$eaf_gwas1,1-dat_merge$eaf_gwas1)
dat_merge$MAF_gwas1[is.na(dat_merge$MAF_gwas1)] <- 0.5
dat_merge$MAF_gwas2 <- ifelse(dat_merge$eaf_gwas2 < 0.5,dat_merge$eaf_gwas2,1-dat_merge$eaf_gwas2)
dat_merge$MAF_gwas2[is.na(dat_merge$MAF_gwas2)] <- 0.5
dat_merge <- dat_merge[dat_merge$MAF_gwas2 != 0, ]


gwas1_form <- list(snp = dat_merge$SNP,
                    position = dat_merge$pos_gwas2,
                    beta = dat_merge$beta_gwas1,
                    varbeta = dat_merge$se_gwas1 ^ 2,
                    MAF = dat_merge$MAF_gwas1,
                    type = "quant",
                    N = dat_merge$N_gwas1)
gwas2_form<- list(snp = dat_merge$SNP,
                   position =dat_merge$pos_gwas2,
                   beta = dat_merge$beta_gwas2,
                   varbeta = dat_merge$se_gwas2^2,
                   MAF=dat_merge$MAF_gwas2,
                   type = "quant",
                   N = dat_merge$N_gwas2)
check_dataset(gwas1_form)
check_dataset(gwas2_form)
my.res <- coloc.abf(dataset1 =gwas1_form ,dataset2 = gwas2_form)

posterior_results <- data.frame(
  nsnps = my.res$summary["nsnps"],
  H0 = my.res$summary["PP.H0.abf"],
  H1 = my.res$summary["PP.H1.abf"],
  H2 = my.res$summary["PP.H2.abf"],
  H3 = my.res$summary["PP.H3.abf"],
  H4 = my.res$summary["PP.H4.abf"]
)
write.csv(posterior_results, file = "coloc_results_eQTL.csv", row.names = FALSE)




#step2 sQTL
#rs10305439
gwas1 <- fread("sqtl_pancreas.txt",header = T)
gwas1 <- gwas1[gwas1$chr == "chr6",]
gwas1 <- gwas1[gwas1$pos> posx-500000   &  gwas1$po < posx+500000,]
gwas1$samplesize<-round(gwas1$mac/gwas1$maf/2)
gwas1<-data.frame(gwas1$snp,
                  gwas1$pos,
                  gwas1$beta,
                  gwas1$se,
                  gwas1$eaf,
                  gwas1$p,
                  gwas1$ref,
                  gwas1$alt,
                  gwas1$phenotype_id,
                  gwas1$samplesize
)
colnames(gwas1)<-c("SNP","pos","beta","se","eaf","pval","other_allele","effect_allele","gene_id","N")
gwas2<-data
gwas2 <- gwas2[gwas2$Chromosome == 6,]
gwas2<- gwas2[gwas2$Position > posy-500000  &  gwas2$Position < posy+500000,]

colnames(gwas2)<-c("SNP","pos","beta","se","eaf","pval","other_allele","effect_allele")
dat_merge <- merge(gwas1,gwas2,by=c("SNP"),suffixes = c("_gwas1","_gwas2"))


dat_merge$label_a1 <- ifelse(dat_merge$effect_allele_gwas1 == dat_merge$effect_allele_gwas2,"yes","no")
dat_merge$label_a1a2 <- ifelse(dat_merge$effect_allele_gwas1 == dat_merge$other_allele_gwas2,"yes","no")
dat_merge$label_a2a1 <- ifelse(dat_merge$other_allele_gwas1 == dat_merge$effect_allele_gwas2,"yes","no")
dat_merge$effect_allele_gwas2 <- ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                         dat_merge$effect_allele_gwas1,
                                         dat_merge$effect_allele_gwas2 )
dat_merge$other_allele_gwas2 <- ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                                        dat_merge$other_allele_gwas1,
                                        dat_merge$other_allele_gwas2 )
dat_merge$beta_gwas2 <-ifelse( dat_merge$label_a1a2 == "yes" & dat_merge$label_a2a1 == "yes",
                               -dat_merge$beta_gwas2 ,
                               dat_merge$beta_gwas2)
dat_merge[,c("label_a1","label_a2","label_a1a2","label_a2a1")] <- NULL

dat_merge$MAF_gwas1   <- ifelse(dat_merge$eaf_gwas1   < 0.5,dat_merge$eaf_gwas1,1-dat_merge$eaf_gwas1)
dat_merge$MAF_gwas1[is.na(dat_merge$MAF_gwas1)] <- 0.5
dat_merge$MAF_gwas2 <- ifelse(dat_merge$eaf_gwas2 < 0.5,dat_merge$eaf_gwas2,1-dat_merge$eaf_gwas2)
dat_merge$MAF_gwas2[is.na(dat_merge$MAF_gwas2)] <- 0.5
dat_merge <- dat_merge[dat_merge$MAF_gwas2 != 0, ]

gwas1_form <- list(snp = dat_merge$SNP,
                   position = dat_merge$pos_gwas2,
                   beta = dat_merge$beta_gwas1,
                   varbeta = dat_merge$se_gwas1 ^ 2,
                   MAF = dat_merge$MAF_gwas1,
                   type = "quant",
                   N = dat_merge$N_gwas1)
gwas2_form<- list(snp = dat_merge$SNP,
                  position =dat_merge$pos_gwas2,
                  beta = dat_merge$beta_gwas2,
                  varbeta = dat_merge$se_gwas2^2,
                  MAF=dat_merge$MAF_gwas2,
                  type = "quant",
                  N = dat_merge$N_gwas2)
check_dataset(gwas1_form)
check_dataset(gwas2_form)
my.res <- coloc.abf(dataset1 =gwas1_form ,dataset2 = gwas2_form)

posterior_results <- data.frame(
  nsnps = my.res$summary["nsnps"],
  H0 = my.res$summary["PP.H0.abf"],
  H1 = my.res$summary["PP.H1.abf"],
  H2 = my.res$summary["PP.H2.abf"],
  H3 = my.res$summary["PP.H3.abf"],
  H4 = my.res$summary["PP.H4.abf"]
)

write.csv(posterior_results, file = "coloc_results_sQTL.csv", row.names = FALSE)













