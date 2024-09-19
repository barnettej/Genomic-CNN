#this code is designed to trim the number of SNPs used in an analysis based on several factors
#here, we select for SNPs with maf > .1, impute_score >= 0.9, and an association p_value of <0.01
#from these SNPs, we then select the SNPs that are most dense with annotations

library(data.table)
library(plyr)


# Load hg19-compatible dataset
setwd("/home/user/user/ukbiobank/t2d/data")
sum_stat <- fread("/home/user/user/ukbiobank/t2d/cid/t2d.sumstat")
colnames(sum_stat)[1] <- "snp"

for(i in 1:22){
  pheno <- "t2d"
  ICD <- "E11"
  run_name <- "t2d.annocut.p0.01"
  
  
  cid <- fread(paste0("/home/user/user/ukbiobank/t2d/data/annotations.c",i,".txt"))
  cid <- cid[,-c(2:6)]
  data = fread(paste0("/home/user/user/ukbiobank/imputation_info/ukb_mfi_chr",i,"_v3.txt"))
  colnames(data) <- c("change","snp","start", "ref","alt","maf","major","impute_score")
  data <- data[data$maf > 0.1 & data$impute_score >= 0.9,]
  cid <- cid[cid$snp %in% data$snp]
  
  # merge sumstat with available snps in cid

  cid <- merge(sum_stat, cid,  by = "snp", sort = F)
  cid <- cid[P<0.01,]
  
  #sum the annotations for each SNP (per row)
  cid$sum <- rowSums(cid[,2:length(cid)], na.rm = T)
  
  
  #look at density plot
  plot(density(cid$sum))
  
  #select rows meeting some threshold of annotations
  threshold <- 2
  cut.cid <- cid[cid$sum > threshold, ]
  
  #plot(density(cut.cid$sum))
  
  #create snp list for trimming data
  snps <- cut.cid[,1]
  write.table(snps, paste0("snps.cid.",run_name,".",i,".txt"), quote = F, col.names = F, row.names = F)
  
  
  #make new bed that only contains the included snps
  bed.loc <- "/mnt/backup/user/gfetch/"
  
  #from scratch
  # system(paste0("plink --bfile ",bed.loc,"ukb_c",i, 
  #               " --extract snps.cid.",run_name,".",i,".txt --make-bed --out trimmed.",run_name,".c",i))
  # system(paste0("plink2 --bfile trimmed.",run_name,".c",i,
  #               " --rm-dup force-first --make-bed --out trimmed.",run_name,".c",i))
  # system(paste0("plink --bfile trimmed.",run_name,".c",i, " --indep-pairwise 50 5 0.8 --out c",i,".",run_name))
  # system(paste0("plink --bfile trimmed.",run_name,".c",i, " --extract c",i,".",run_name,
  #               ".prune.in --make-bed --out trimmed.",run_name,".c",i))
  # 
  #or from larger p
  system(paste0("plink --bfile trimmed.c",i, " --extract snps.cid.",
                run_name,".",i,".txt --indep-pairwise 50 5 0.8 --make-bed --out trimmed.",run_name,".c",i))
}

  #create snp list of trimmed data
  system(paste0("cat trimmed.",run_name,".c*.prune.in > snps.cid.",run_name,".txt"))
  

  