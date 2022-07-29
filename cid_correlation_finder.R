
library(data.table)
library(plyr)


# Load hg19-compatible dataset
setwd("/home/user/eric/ukbiobank/t1d/data2")


for(i in 22:1){
  start_time <- Sys.time()
  merge.revel <- fread(paste0("/home/user/eric/ukbiobank/t1d/data/annotations.c",i,".txt"))
  bed.loc <- "/home/user/eric/ukbiobank/t1d/data/"
  
  #### highest correlation with annotation #### 
  ## determine whether each annotation is binary or continuous
  annot.info <- merge.revel[,stack(colSums(merge.revel[-1] > 1) == 0)]
  binary.annot <- merge.revel[,which(stack(colSums(merge.revel[-1] > 1) == 0) == T)]
  
  
  
  
  
  r2.annots <- as.data.table(merge.revel[,snp])
  colnames(r2.annots) <- "snp"
  for( k in binary.annot) {
    include <- c("snp",colnames(merge.revel)[k])
    annot <- merge.revel[, ..include]
    backup.colnames <- colnames(annot)
    colnames(annot) <- c("snp", "annot")
    snps.with.annot <- annot[annot == 1, "snp"]
    
    annot.r2 <- data.table()
    for(j in 1:dim(merge.revel)[1]){
      snp <- as.data.table(annot$snp[j])
      if(annot$annot[j] == 1){
        annot.r2.j <- data.table("SNP" = snp, "annot.r2" = 1)
      } else {
        snp.list <- rbind(snp, snps.with.annot$snp, use.names = F)
        write.table(snp.list, "snp_list.txt", col.names = F, row.names = F, quote = F)
        system(paste0("plink  --memory 10000 --bfile ",bed.loc,"c",i,".10k  --extract snp_list.txt --r2 --ld-window 20 --ld-window-kb 1000 --ld-window-r2 0 --ld-snp ", snp))
        r2 <- fread("plink.ld")
        if(dim(r2)[1] > 0) {
          r2.with.annot <- merge.data.table(snps.with.annot, r2, by.x = "snp", by.y = "SNP_B")
        } else {
          r2.with.annot <- data.table()
        }
        if(dim(r2.with.annot)[1] > 0) {
          annot.r2.j <- data.table("SNP" = snp, "annot.r2" = max(r2.with.annot$R2))
        } else { #this next line adds a 0 if there are no correlated snps with the annotation
          annot.r2.j <- data.table("SNP" = snp, "annot.r2" = 0)
        }}
      annot.r2 <- rbind(annot.r2, annot.r2.j, use.names = F)
    }
    colnames(annot.r2) <- backup.colnames
    colnames(annot.r2)[2] <- paste0(colnames(merge.revel)[k],".r2")
    r2.annots <- merge.data.table(r2.annots, annot.r2)
    r2.annots <- r2.annots[!duplicated(r2.annots$snp),]
  }
  
  
  write.table(r2.annots, paste0("annotations.c",i,".r2.txt"), col.names = T, row.names = F)
  
  cid <- merge.data.table(merge.revel, r2.annots, by = "snp")
  
  
  
  ##### fill in missing values with value of most correlated SNP  
  
  cont.annot <- merge.revel[,which(is.na(stack(colSums(merge.revel[-1] > 1) == 0) == T))]
  
  combined.cont.annot <- as.data.table(merge.revel[,snp])
  colnames(combined.cont.annot) <- "snp"
  for( j in cont.annot) {
    include <- c("snp",colnames(merge.revel)[j])
    annot <- merge.revel[, ..include]
    backup.colnames <- colnames(annot)
    colnames(annot) <- c("snp", "annot")
    snps.with.annot <- annot[complete.cases(annot),]
    na.annot <- annot[!complete.cases(annot),]
    annot.replace <- na.annot
    annot.replace$replacement.snp <- ""
    for(k in 1:dim(na.annot)[1]){
      snp <- na.annot$snp[k]
      snp.list <- rbind(snp, snps.with.annot$snp, use.names = F)
      write.table(snp.list, "snp_list.txt", col.names = F, row.names = F, quote = F)
      system(paste0("plink --memory 10000 --bfile ",bed.loc,"c",i,
                    ".10k --r2 --extract snp_list.txt --ld-window 20 --ld-window-kb 10000 --ld-window-r2 0 --ld-snp ",
                    snp, " --out na_replace"))
      r2 <- read.table("na_replace.ld", header = TRUE)
      r2.with.annot <- merge(snps.with.annot, r2, by.x = "snp", by.y = "SNP_B")
      #annot.r2 <- data.table("SNP" = snp, "annot.r2" = max(r2$V2))
      annot.replace$replacement.snp[k] <- r2.with.annot$snp[which.max(r2.with.annot$R2)]
      
    }
    #annot.replace <- merge(annot.replace, snps.with.annot, by.x = "replacement.snp", by.y = "snp")
    annot.replace <- merge(annot.replace, snps.with.annot, by.x = "replacement.snp", by.y = "snp")
    
    annot <- annot[is.na(annot), annot := annot.replace[.SD, annot.y, on="snp"]]
    colnames(annot) <- backup.colnames
    colnames(annot)[2] <- paste0(colnames(merge.revel)[j],".na_replaced")
    combined.cont.annot <- merge.data.table(combined.cont.annot, annot)
    gc()
  }
  
  
  
  write.table(combined.cont.annot, paste0("annotations.c",i,".na_rep.txt"), col.names = T, row.names = F)
  
  cid <- merge.data.table(cid, combined.cont.annot, by = "snp")
  
  write.table(cid, paste0("cid.c",i,".txt"), col.names = T, row.names = F)
  
  end_time <- Sys.time()
  time_taken <- (end_time - start_time)
  print(time_taken)
  
  writeLines(paste("chr",i,"took",time_taken), paste0("chr",i,".time"))
  rm(list = ls())
  gc()
  
}


