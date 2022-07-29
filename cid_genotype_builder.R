
library(data.table)
library(plyr)
library(tidyr)
pheno <- "t1d"
ICD <- "E10"
run_name <- "t1d.sumstat.p0.1"

setwd("/home/user/eric/ukbiobank/t1d/cid")
full.cid <- data.table()
full.cid <- fread("cid.c1.txt")
full.cid <- full.cid[,-21]
for( i in 2:22){
  chr.cid <- fread(paste0("cid.c",i,".txt"))
  common.cols <- intersect(colnames(full.cid), colnames(chr.cid))
  full.cid <- rbind(
    subset(full.cid, select = common.cols),
    subset(chr.cid, select = common.cols))
}
  
sum_stat <- fread("/home/user/eric/gcnn/cid_files/converted_t1d_sumstat")
colnames(sum_stat)[1] <- "snp"

#### scale data between 0 and 1
fn <- function(x) (1 - 0)*((x-min(x))/(max(x)-min(x)))+ 0
#this is to include all annotations
#cid <- data.frame(lapply(full.cid[,7:length(full.cid)], fn))
#this only includes the r2 and na.replaced features
cid <- data.frame(lapply(full.cid[,18:length(full.cid)], fn))
cid <- cbind(full.cid[,1],cid)
cid <- merge(sum_stat[,c(1,6,8)], cid,  by = "snp")

#remove revel for now
#cid <- cid[,-14]

cut.cid <- cid[P<0.1,]
write.table(cut.cid, paste0(run_name,".scaled.cid.annotations.txt"), quote = F, col.names = T, row.names = F)
# # #look at density plot
# plot(density(cid$sum))
# 
# #save plot
# jpeg(file="annotation.density_plot.jpeg")
# plot(density(cid$sum))
# dev.off()
# 
# #select rows meeting some threshold of annotations
# threshold <- 4
# cut.cid <- cid[cid$sum > threshold, ]
# 
# plot(density(cut.cid$sum))

#create snp list for trimming data
snps <- cut.cid[,1]
write.table(snps, paste0("snps.cid.",run_name,".txt"), quote = F, col.names = F, row.names = F)

#make new bed that only contains the included snps
bed.loc <- "/mnt/backup/user/gfetch/"
for(i in 1:22){
  system(paste0("plink --bfile ",bed.loc,"ukb_c",i, 
                " --extract snps.cid.",run_name,".txt --make-bed --out trimmed.c",i))
  system(paste0("plink2 --bfile trimmed.c",i,
                " --rm-dup force-first --make-bed --out trimmed.c",i))
}


#combine all chr and extract snps

system(paste0("plink --merge-list merge_file.txt --make-bed ",
              "--extract snps.cid.",run_name,".txt",
              " --out merged"))
system(paste0("plink --bed merged.bed --bim merged.bim --fam train_split.fam",
              " --must-have-sex --make-bed --out cid"))

## phenotype file generated from ukbiobank_phenotype.R

## matched pairing
subsets <- c("train","valid","test")
for(i in subsets){
  #check for imbalances
  df <- fread("t1d.covar.txt")
  df <- df[df$V1 == i,]
  colnames(df)[6] <- "t1d"
  df$t1d[df$t1d == 1] <- 0
  df$t1d[df$t1d == 2] <- 1
  m.out0 <- matchit(t1d ~ scaled_age + sex, data = df, method = NULL, distance = "glm")
  summary(m.out0)
  #match
  m.out1 <- matchit(t1d ~ scaled_age + sex, data = df, method = "nearest", distance = "glm")
  #test balance of match
  summary(m.out1, un = F)
  #get output
  matched <- match.data(m.out1)
  write.table(matched, paste0("t1d.covar.",i,".matched.txt"), quote = F, col.names = T, row.names = F)
  matched$t1d[matched$t1d == 1] <- 2
  matched$t1d[matched$t1d == 0] <- 1
  matched <- setcolorder(matched, c("V1", "V2"))
  write.table(matched, paste0("t1d.covar.",i,".matched.fam"), quote = F, col.names = F, row.names = F)
}
system(paste0("cat t1d.covar.*.matched.fam > t1d.covar.matched.fam"))


system(paste0("plink --bed cid.bed --bim cid.bim --fam cid.fam --keep ",pheno,".covar.matched.fam",
              " --make-bed --out matched"))

#### make subset input for p value trimmed version####
subsets <- c("train","valid","test")
for(i in subsets){

  
  system(paste0("plink --bed matched.bed --bim matched.bim --fam ",pheno,".covar.matched.fam --extract snps.cid.",run_name,".txt --prune --must-have-sex --family ",
                "-keep-cluster-names ",i," --make-bed --out ",run_name,".",i))
  
  # recode to file format that includes allele count for each person
  system(paste0("plink --bfile ",run_name,".",i," --recode structure --family -keep-cluster-names ",
                i," --out ",run_name,".",i))
  # split file into smaller bits for next part
  system(paste0("split -l 10000 -d ",run_name,".",i,".recode.strct_in ",i,".segment."))
  
  #load in phenotype fam (made in ukbiobank_phenotype.R) and change for tf formatting
  fam <- fread(paste0(run_name,".",i,".fam"))
  fam$V6[fam$V6 == 1] <- 0
  fam$V6[fam$V6 == 2] <- 1
  colnames(fam)[6] <- pheno
  
  #merge fam with age and sex info
  age.sex <- fread("age+sex.csv")
  fam <- merge(fam, age.sex, by.x = "V2", by.y = "f.eid")
  
  
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
      plink.recode <- cbind(fam[c(1:max.people),c(6,7,8)], plink.recode)
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
      plink.recode <- cbind(fam[c(min.people:max.people),c(6,7,8)], plink.recode)
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
