library(MatchIt)
library(data.table)
library(plyr)
library(tidyr)
pheno <- "t2d"
run_name <- "t2d"

setwd("/home/user/user/ukbiobank/t2d/data")


#create snp list of trimmed snps and snps that were pruned due to ld
system(paste0("cat trimmed.",run_name,".c*.prune.out > snps.pruned.",run_name,".txt"))
system(paste0("cat snps.pruned.",run_name,".txt snps.cid.",
              run_name,".txt > snps.within_ld.",run_name,".txt"))
system(paste0("cat t2d.match_list.*.txt > t2d.match_list.full.txt"))

for(k in 1:22){
    bed.loc <- "/mnt/backup/user/gfetch/"
    system(paste0("plink --make-bed --maf 0.3 --snps-only just-acgt --bfile ",bed.loc,"ukb_c",k,
                  " --pheno t2d.pheno --exclude snps.within_ld.",
                  run_name,".txt --keep t2d.match_list.full.txt",
                  " --out ancestry.c",k))

}
system(paste0("plink --merge-list ancestry_merge_file.txt  --make-bed ",
                " --out ",run_name,".ancestry.full"))
  #if this results in multiallelic error make new exclude file
system(paste0("cat snps.within_ld.",run_name,".txt ",
                run_name,".ancestry.full-merge.missnp > combined_exclude.txt"))

for(k in 1:22){
    bed.loc <- "/mnt/backup/user/gfetch/"
    system(paste0("plink --make-bed --maf 0.3 --bfile ancestry.c",k,
                  " --pheno t2d.pheno --set-missing-var-ids @:#[b37]\\$1,\\$2 --exclude combined_exclude.txt",
                  " --keep t2d.match_list.full.txt",
                  " --out ancestry.c",k))
    
}
  
system(paste0("plink --merge-list ancestry_merge_file.txt  --make-bed ",
                " --out ",run_name,".ancestry.full"))
  

#run pca 
system(paste0("plink2 --memory 10000 --bed ",run_name,".ancestry.full.bed --bim ",
              run_name,".ancestry.full.bim --fam train_split.fam --indep-pairwise 100 10 0.2 --maf .10 --out plink_files/ancestry.pca"))
system(paste0("plink --memory 10000 --bed ",run_name,".ancestry.full.bed --bim ",
              run_name,".ancestry.full.bim --fam train_split.fam  --extract plink_files/ancestry.pca.prune.in --family --pca-cluster-names train ",
              "--pca  --out plink_files/ancestry.pca"))

pca <- fread("plink_files/ancestry.10.pca.eigenvec")
pca.s <- as.data.table(scale(pca[,c(3:12)]))
pca <-cbind(pca[,1],pca.s)
colnames(pca)[1] <- "f.eid"
subsets <- c("train","valid","test")
for(i in subsets){
  input <- fread(paste0(run_name,".",i,".input.csv"))
  merged <- merge(pca, input[,-c(3,4)], by = "f.eid")
  merged <- merged[sample(nrow(merged)),]
  write.table(merged[,-1], paste0(run_name,".",i,".input.with_pca.scale.csv"),sep = ",", quote = F, row.names = F, col.names = T)
  write.table(merged[,c(2:12)], paste0(run_name,".",i,".pca.scale.csv"),sep = ",", quote = F, row.names = F, col.names = T)
}

df <- fread(paste0(run_name,".test.pca.csv"))
df.scale <- as.data.table(scale(df))