## MR analysis testing causal effects of GLP-1R effect on outcomes selected
rm(list=ls())
gc()
library(TwoSampleMR)
library(coloc)
library(vroom)


## Step 1:  Constructing exposure data with instruments validation 
data<-read.csv("multi-omic_qtl_instruments.csv")
#exposure data
data <- data.frame(SNP = data$SNP,
                   chr = data$chr,
                   pos = data$pos,
                   beta = as.numeric(data$beta),
                   se = as.numeric(data$se),
                   effect_allele = data$effect_allele,
                   other_allele = data$other_allele,
                   eaf=data$effect_allele_freq,
                   pval = as.numeric(data$pval),
                   Phenotype = data$Phenotype,
                   samplesize = data$samplesize,
                   info = paste(data$Study,data$cis_trans,sep="_"))
data <- format_data(data, type="exposure", phenotype_col = "Phenotype",chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize",info_col = "info")
#ld clumping
data <-clump_data(data,clump_kb=10000,clump_r2=0.001,clump_p1=1,clump_p2=1,pop="EUR")
# F statistics>10 
# Calculate F-statistics using R2 and MAF
exposure_dat <- cbind(data, fstatistics = 1)
for (s in 1:nrow(exposure_dat)){
  beta <- exposure_dat[s, "beta.exposure"]
  n <- exposure_dat[s, "samplesize.exposure"]
  eaf <- exposure_dat[s, "eaf.exposure"]
  
  # Calculate R2 using the formula: R2 = beta^2 * 2 * MAF * (1-MAF)
  r2 <- beta^2 * 2 * eaf * (1 - eaf)
  
  # Calculate F-statistic
  exposure_dat[s, "fstatistics"] <- (r2 * (n - 2)) / (1 - r2)
}

print(min(exposure_dat$fstatistics))
print(max(exposure_dat$fstatistics))
exposure_dat <- exposure_dat[exposure_dat$fstatistics > 10, ]


## Step 2:  MR analysis testing causal effects of GLP-1R expression on positive control/overall ovarian cancer/EOC subtypes
#outcome data 
outcomelist<-read.csv("outcomelist.csv")
i<-1
out_data<- vroom(outcomelist[i,"filename"])
outcome_dat<-format_data(
    out_data,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se", 
    eaf_col = "eaf",
    effect_allele_col ="effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col= "samplesize",
    chr_col='chr',
    pos_col='pos'
  )
#harmonization
har_dat <- harmonise_data(exposure_dat, outcome_dat, action=2)
#steiger filtering
stei_dat <- steiger_filtering(har_dat)
#MR results  
mr_results <- mr(har_dat)
mr_heterogeneity(har_dat)
mr_pleiotropy_test(har_dat)

#rerun Step 2 after replacing i from 2 to nrow(outcomelist)








