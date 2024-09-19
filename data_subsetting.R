library(data.table)
library(plyr)
library(tidyr)
#make fam for data splitting
setwd("/home/user/user/ukbiobank/data")

system("plink --bfile ukb_c22 --make-just-fam --out train_split")


## split data 
train_split <- 0.70
valid_split <- 0.15
test_split <- 0.15

#read in fam for lengths
df <- fread("train_split.fam", header = F)
train_split <- train_split*nrow(df)
valid_split <- valid_split*nrow(df)

#use this file if FID isn't duplicate of IID
#df$V2 <- paste0(df$V1,"+",df$V2)
df$V1 <- sample(1:nrow(df), size = nrow(df), replace = F)
df$V1[df$V1 > (train_split + valid_split)] <- "test"
df$V1[as.numeric(df$V1) > train_split & as.numeric(df$V1) <= (train_split + valid_split)] <- "valid"
df$V1[as.numeric(df$V1) <= train_split] <- "train"

write.table(df, "train_split.fam", col.names = F, row.names = F, quote = F)




#make 10k subset fam file for correlation analyses 
i = 22
system(paste0("plink --bed ukb_c",i,".bed --bim ukb_c",i,".bim --fam train_split.fam --make-just-fam",
              " --family --keep-cluster-names train --out c",i,".train"))
train <- fread(paste0("c",i,".train.fam"))
train <- train[sample(nrow(train), 10000),]
write.table(train, "train.10k.fam", quote = F, col.names = F, row.names = F)


