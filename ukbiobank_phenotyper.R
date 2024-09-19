library(tidyr)
library(data.table)
library(MatchIt)
setwd("/home/user/user/ukbiobank/t2d/data")
bd <- fread("/home/user/user/ukbiobank/gfetch/ukb.tab")
pheno.list <- grep(pattern = "f.41270.0.", x = colnames(bd))
pheno.list <- c(1, pheno.list)
pheno <- subset(bd, select = pheno.list)
rm(bd)
gc()
eid <- subset (pheno, select = 1)
pheno.c <- pheno %>%
  unite('Merged', f.41270.0.0:f.41270.0.242, ,sep = ",", remove = TRUE)

have.pheno <- pheno.c[pheno.c$Merged %like% "E11", ]
have.pheno$pheno <- 1
pheno <- merge(pheno.c, have.pheno, by = c("f.eid", "Merged"), all.x = TRUE, sort = F)
pheno$pheno[is.na(pheno$pheno)] <- 0
pheno <- pheno[,-2]
pheno.file <- data.table("fid" = pheno$f.eid, "iid" = pheno$f.eid, "t2d" = pheno$pheno)
pheno.file$t2d <- pheno.file$t2d + 1
write.table(pheno.file, "t2d.pheno", quote = F, col.names = F, row.names = F)

age.sex <- fread("../cid/age+sex.csv")
merge <- merge(pheno, age.sex, by = "f.eid")
colnames(merge)[2] <- "t2d"
write.table(merge, "t2d.covar.txt", quote = F, col.names = T, row.names = F)


## matched pairing
subsets <- c("train","valid","test")
for(i in subsets){
  train.split <- fread("train_split.fam")
  train.split <- as.data.table(train.split$V2[train.split$V1 == i])
  #check for imbalances
  df <- fread(paste0("t2d.covar.txt"))
  df <- merge(df, train.split, by.x = "f.eid", by.y = "V1", sort = F)
  m.out0 <- matchit(t2d ~ scaled_age + sex, data = df, method = NULL, distance = "glm")
  summary(m.out0)
  #match (ratio means we're doing 1:1 matching)
  m.out1 <- matchit(t2d ~ scaled_age + sex, data = df, method = "nearest", distance = "glm", ratio = 1)
  #test balance of match
  summary(m.out1, un = F)
  #get output
  matched <- match.data(m.out1)
  output <- data.table("fid" = matched$f.eid, "iid" = matched$f.eid)
  write.table(output, paste0("t2d.match_list.",i,".txt"), quote = F, col.names = F, row.names = F)
  
  unmatched <- match.data(m.out1, drop.unmatched = F)
  unmatched <- unmatched[unmatched$weights==0]
  unmatched <- data.table("fid" = unmatched$f.eid, "iid" = unmatched$f.eid)
  write.table(unmatched, paste0("t2d.unmatch_list.",i,".txt"), quote = F, col.names = F, row.names = F)
  
}


system(paste0("cat t2d.match_list.*.txt > t2d.match_list.full.txt"))
system(paste0("cat t2d.unmatch_list.*.txt > t2d.unmatch_list.full.txt"))


