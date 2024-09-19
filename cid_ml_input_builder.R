library(MatchIt)
library(data.table)
library(plyr)
library(tidyr)
pheno <- "t2d"
ICD <- "E11"
run_name <- "t2d"
icd.list <- c("E66", "E78", "I10", "I25", "K80", "I20", "N39", "K29", "R07")
setwd("/home/user/user/ukbiobank/t2d/data")

full.cid <- data.table()
full.cid <- fread("cid.c1.txt")
data = fread(paste0("/home/user/user/ukbiobank/imputation_info/ukb_mfi_chr1_v3.txt"))
colnames(data) <- c("change","snp","start", "ref","alt","maf","major","impute_score")
data <- data[data$maf > 0.1 & data$impute_score >= 0.9,]
full.cid <- full.cid[full.cid$snp %in% data$snp]

for( i in 2:22){
  chr.cid <- fread(paste0("cid.c",i,".txt"))
  data = fread(paste0("/home/user/user/ukbiobank/imputation_info/ukb_mfi_chr",i,"_v3.txt"))
  colnames(data) <- c("change","snp","start", "ref","alt","maf","major","impute_score")
  data <- data[data$maf > 0.1 & data$impute_score >= 0.9,]
  chr.cid <- chr.cid[chr.cid$snp %in% data$snp]
  common.cols <- intersect(colnames(full.cid), colnames(chr.cid))
  full.cid <- rbind(
    subset(full.cid, select = common.cols),
    subset(chr.cid, select = common.cols))
}


sum_stat <- fread("/home/user/user/ukbiobank/t2d/cid/t2d.sumstat")
colnames(sum_stat)[1] <- "snp"

#### scale data between 0 and 1
fn <- function(x) (1 - 0)*((x-min(x))/(max(x)-min(x)))+ 0
#this is to include all annotations
#cid <- data.frame(lapply(full.cid[,7:length(full.cid)], fn))
#this only includes the r2 and na.replaced features
cid <- data.frame(lapply(full.cid[,34:length(full.cid)], fn))
cid <- cbind(full.cid[,1],cid)
cid <- merge(sum_stat, cid,  by = "snp", sort = F)

#read in snp list if used cid_snp_trimmer
cut.snps <- fread(paste0("snps.cid.",run_name,".txt"), header = F)
cut.cid <- merge(cid, cut.snps, by.x = "snp", by.y = "V1", sort = F)


# if you didn't use cid_snp_trimmer, create snp list of trimmed data
# cut.cid <- cid[P<0.1,]
# system(paste0("cat c*.prune.in > snps.cid.",run_name,".txt"))


## attach correlated phenotypes
setwd("/home/user/user/ukbiobank/t2d/data")
for(k in 1:22){
  bed.loc <- "/mnt/backup/user/gfetch/"
  system(paste0("plink --make-bed --bfile ",bed.loc,"ukb_c",k,
                " --pheno t2d.pheno --extract trimmed.",run_name,".c",k,".prune.in",
                " --keep t2d.unmatch_list.full.txt",
                " --out trimmed.unmatched.c",k))
}
system(paste0("plink --merge-list merge_file.unmatched.txt  --make-bed ",
              " --out ",run_name,".unmatched"))
for(k in 1:22){
  bed.loc <- "/mnt/backup/user/gfetch/"
  system(paste0("plink --make-bed --bfile trimmed.unmatched.c",k," --exclude ",
                run_name,".unmatched-merge.missnp --extract trimmed.",run_name,".c",k,
                ".prune.in --keep t2d.unmatch_list.full.txt",
                " --out trimmed.unmatched.c",k))
}
system(paste0("plink --merge-list merge_file.unmatched.txt  --make-bed ",
              " --out ",run_name,".unmatched"))


unmatched.fam <- fread(paste0(run_name,".unmatched.fam"))

for(i in icd.list){
  pheno.df <- fread("/home/user/user/ukbiobank/pheno.df.txt")
  #attach phenotype
  eid <- subset (pheno.df, select = 1)
  pheno.c <- pheno.df %>%
    unite('Merged', 2:length(pheno.df),sep = ",", remove = TRUE)
  
  have.pheno <- pheno.c[pheno.c$Merged %like% i, ]
  if(dim(have.pheno)[1]>=1000){
    have.pheno$pheno <- 1
    pheno.df <- merge(pheno.c, have.pheno, by = c("f.eid", "Merged"), all.x = T, sort = F)
    pheno.df$pheno[is.na(pheno.df$pheno)] <- 0
    pheno.df <- pheno.df[,-2]
    
    new.fam <- merge(unmatched.fam, pheno.df, by.x = "V2", by.y = "f.eid", all.x = T, sort = F)
    new.fam$pheno <- new.fam$pheno + 1
    pheno.file <- data.table("fid" = new.fam$V2, "iid" = new.fam$V2, "pheno" = new.fam$pheno)
    write.table(pheno.file, paste0(i,".unmatched.pheno"), quote = F, col.names = F, row.names = F)
  }
  rm(new.fam, have.pheno, eid, pheno.c, pheno.file,pheno.df)
  gc()
  
}

icd.list <- list.files(path = ".", pattern = "*.unmatched.pheno")
assoc.loc <- "/home/user/user/ukbiobank/t2d/cid/"
for(i in icd.list){
  ##gwas with just unmatched
  system(paste0("plink --bfile ",run_name,".unmatched --pheno ",i," --assoc --out ",assoc.loc,i))
  
}


setwd("/home/user/user/ukbiobank/t2d/data")

for(i in icd.list){
  risk <- fread(paste0(assoc.loc,i,".assoc"))
  colnames(risk)[2] <- "snp"
  colnames(risk)[10] <- "OR"
  risk <- risk[,c(2,10)]
  risk$risk <- log10(risk$OR)
  risk <- risk[,c(1,3)]
  colnames(risk)[2] <- i
  cut.cid <- merge.data.table(risk, cut.cid, by = "snp", sort = F)
}

write.table(cut.cid, paste0(run_name,".scaled.cid.annotations.txt"), quote = F, col.names = T, row.names = F)


## phenotype file generated from ukbiobank_phenotyper.R


#### make subset input for p value trimmed version####

subsets <- c("train","valid","test")
for(i in subsets){
  for(k in 1:22){
    bed.loc <- "/mnt/backup/user/gfetch/"
    system(paste0("plink --make-bed --bfile ",bed.loc,"ukb_c",k," --pheno t2d.pheno --exclude ",
                  run_name,".unmatched-merge.missnp --extract snps.cid.",run_name,".txt --keep t2d.match_list.",i,".txt",
                  " --out trimmed.c",k))
  }
  system(paste0("plink --merge-list merge_file.txt  --make-bed ",
                " --out ",run_name,".",i))
  #if this results in multiallelic error make new exclude file

  # recode to file format that includes allele count for each person
  system(paste0("plink --bfile ",run_name,".",i," --recode structure --keep t2d.match_list.",i,".txt",
                " --out ",run_name,".",i))
  # split file into smaller bits for next part
  system(paste0("split -l 10000 -d ",run_name,".",i,".recode.strct_in ",i,".segment."))
  
  #load in phenotype fam (made in ukbiobank_phenotype.R) and change for tf formatting
  fam <- fread(paste0(run_name,".",i,".fam"))
  fam$V6[fam$V6 == 1] <- 0
  fam$V6[fam$V6 == 2] <- 1
  colnames(fam)[6] <- pheno
  colnames(fam)[2] <- "f.eid"
  #merge fam with age and sex info
  age.sex <- fread(paste0(assoc.loc,"age+sex.csv"))
  fam <- merge(fam, age.sex, by.x = "f.eid", by.y = "f.eid", sort = F)
  
  
  #loop through all segments of structure file and output/append input for tf
  structure.list <- list.files(path="./", pattern = paste0(i,".segment.[0-9][0-9]$"))
  for(k in structure.list){
    if(k == structure.list[1]){
      ## load in snp id and combine with summed alleles
      plink.recode.snps <- read.table(paste0(i,".segment.00"), nrows = 1)
      plink.recode <- fread(k, skip = 2)
      #get list of IDs
      plink.recode.id <- plink.recode[,c(1,2)]
      ## change 2 (representing A2 = major allele) to 0 so we can get a count of minor alleles
      plink.recode[,3:length(plink.recode)][plink.recode[,3:length(plink.recode)] == 2] <- 0
      plink.recode <- as.matrix(plink.recode[,3:length(plink.recode)])
      #plink.recode <- as.matrix(plink.recode[,3:length(plink.recode)]-1)
      plink.recode <- as.data.table(plink.recode[, seq(1, ncol(plink.recode),
                                                       by = 2)] + plink.recode[, seq(2, ncol(plink.recode), by = 2)])
      #need to keep track of number of people in each segment so that we cbind to correct label
      max.people <- dim(plink.recode)[1]
      gc()
      colnames(plink.recode) <- as.character(plink.recode.snps)
      plink.recode <- cbind(fam[c(1:max.people),c(1,6,7,8)], plink.recode)
      plink.recode.snps <- transpose(plink.recode.snps)
      colnames(plink.recode.snps) <- "snp"
      write.table(plink.recode, paste0(k,".input.",run_name,".csv"), sep = ",", quote = F, col.names = T, row.names = F)
      #merge with annotations to make sure snps are in same position in both parts
      merged.annotations <- fread(paste0(run_name,".scaled.cid.annotations.txt"))
      merged.annotations <- merged.annotations[!duplicated(merged.annotations$snp),]
      merged.alleles <- merge(plink.recode.snps, merged.annotations, sort = F, by = "snp")
      merged.annotations <- merged.alleles[,c(1,(length(merged.alleles)-(length(merged.annotations)-2)):(length(merged.alleles)))]
      write.table(merged.annotations, paste0(run_name,".",i,".annotations.txt"), col.names = T, row.names = F, quote = F)
      
      transposed.annotations <- transpose(merged.annotations)
      write.table(transposed.annotations, paste0(run_name, ".",i,".input.cid.annotations.txt"), col.names = F, row.names = F, quote = F)
      write.table(transposed.annotations, paste0(run_name, ".",i,".input.cid.annotations.csv"), sep = ",", row.names = F, quote = F, col.names = F)
      
      
      rm(plink.recode)
      gc()
    } else {
      #new minimum should be 1 above the previous maximum
      min.people = max.people + 1
      #load in other segments
      plink.recode <- fread(k)
      #new maximum is the new minimum + the number of people in this segment
      max.people = max.people + dim(plink.recode)[1]
      ## change 2 (representing A2 = major allele) to 0 so we can get a count of minor alleles
      plink.recode[,3:length(plink.recode)][plink.recode[,3:length(plink.recode)] == 2] <- 0
      plink.recode <- as.matrix(plink.recode[,3:length(plink.recode)])
      #plink.recode <- as.matrix(plink.recode[,3:length(plink.recode)]-1)
      plink.recode <- as.data.table(plink.recode[, seq(1, ncol(plink.recode),
                                                       by = 2)] + plink.recode[, seq(2, ncol(plink.recode), by = 2)])
      gc()
      plink.recode <- cbind(fam[c(min.people:max.people),c(1,6,7,8)], plink.recode)
      write.table(plink.recode, paste0(k,".input.",run_name,".csv"), sep = ",", quote = F, col.names = F, row.names = F)
      rm(plink.recode)
      
      gc()
    }
    print(paste0(k, " done."))
  }
  system(paste0("cat ",i,".segment.*.input.",run_name,".csv > ",run_name,".",i,".input.csv"))
}



system(paste0("rm ",run_name,".",i,".recode.strct_in"))
system("rm *segment*")
