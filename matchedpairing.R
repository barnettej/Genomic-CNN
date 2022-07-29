library(data.table)
library(caret)
library(MatchIt)

setwd("/home/user/eric/ukbiobank/cid2/cid")
run_name <- "t1d.sumstat.0.5"


## matched pairing
subsets <- c("train","valid","test")
unmatched.all <- data.table()
for(i in subsets){
#check for imbalances
  df <- fread(paste0(run_name,".",i,".input.csv"))
  m.out0 <- matchit(t1d ~ scaled_age + sex, data = df, method = NULL, distance = "glm")
  summary(m.out0)
  #match
  m.out1 <- matchit(t1d ~ scaled_age + sex, data = df, method = "nearest", distance = "glm")
  #test balance of match
  summary(m.out1, un = F)
  #get output
  matched <- match.data(m.out1)
  write.table(matched, paste0(run_name,".",i,".input.matched.csv"), sep = ",", quote = F, col.names = T, row.names = F)
  ## combine unmatched
  unmatched = df[!row.names(df) %in% row.names(matched),]
  unmatched.all <- rbind(unmatched.all, unmatched)
  }
