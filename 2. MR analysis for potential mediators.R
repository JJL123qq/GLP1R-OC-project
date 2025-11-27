##MR analysis for potential mediators


## Step 1: MR analysis testing causal effects of GLP1R bioactivity on BMI/sex hormone
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

#outcome data
outcomelist<-read.csv("outcomelist_mediators.csv")
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
#run the above codes with i from 1 to 7
#remove the instruments of GLP1R bioactivity proxied by pQTLs when the mediator is protein hormone
#save the above results


## Step 2: MR analysis testing causal effects of  BMI/sex hormone  on subtypes of EOC
exposurelist<-read.csv("outcomelist_mediators.csv")
i<-1
exp_data<- vroom(outcomelist[i,"filename"])
exposure_dat<-format_data(
  exp_data,
  type = "exposure",
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
outcomelist<-c("ieu-a-1123","ieu-a-1124","ieu-a-1125","ieu-a-1121","ieu-a-1122")
har_all<data.frame()
mr_all<-data.frame()
stei_all<-data.frame()

for (m in c(1:5)) {
  outcome_dat<- extract_outcome_data(
    snps = exposure_dat$SNP,
    outcomes = outcomelist[m])
  har_dat <- harmonise_data(
    exposure_dat, 
    outcome_dat
  ) 
  #steiger filtering
  stei_dat <- steiger_filtering(har_dat)
  #MR results
  mr_dat <- mr(har_dat)
  het<-mr_heterogeneity(har_dat)
  print(het)
  ple<-mr_pleiotropy_test(har_dat)
  print(ple)
  
  har_all<-rbind(har_all,har_dat)
  mr_all<-rbind(mr_all,mr_dat)
  stei_all<-rbind(stei_all,stei_dat)
  
}

#run the above codes with i from 1 to 7
#save the above results
#count the weight of mediators in the association of GLP1R and subtypes of EOC

