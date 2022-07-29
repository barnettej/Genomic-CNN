
#library(readr)
library(AnnotationHub)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(plyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)

ah <- AnnotationHub()
grs <- ah[ah$rdataclass == "GRanges" & ah$genome == "hg19",]


# Load hg19-compatible dataset
setwd("/home/user/eric/ukbiobank/t1d/data")
bed.loc <- "/mnt/backup/user/gfetch/"

#txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene



for(i in 21:1){
  ah <- AnnotationHub()
  grs <- ah[ah$rdataclass == "GRanges" & ah$genome == "hg19",]
  
  
  # Load hg19-compatible dataset
  setwd("/home/user/eric/ukbiobank/t1d/data")
  bed.loc <- "/mnt/backup/user/gfetch/"
  
  #txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  
  start_time <- Sys.time()
  data = fread(paste0("/home/user/eric/ukbiobank/imputation_info/ukb_mfi_chr",i,"_v3.txt"))
  bim = fread(paste0(bed.loc,"ukb_c",i,".bim"))
  colnames(bim) <- c("chr","snp","cm","start","alt","ref")
  colnames(data) <- c("change","snp","start", "ref","alt","maf","major","impute_score")
  
  
  data <- data[data$maf > 0.05 & data$impute_score >= 0.9,]
  #load in summary statistics file so we can match with our data
  sum_stat <- fread("/home/user/eric/gcnn/cid_files/converted_t1d_sumstat")
  colnames(sum_stat)[1] <- "snp"
  sum_stat <- sum_stat[,1]
  #merge and keep only snps that are found in both
  data <- merge(data, sum_stat, by = "snp")
  snps <- data[,1]
  write.table(snps, paste0("snps.c",i), quote = F, col.names = F, row.names = F)
  
  ## subset dataset for correlation finding
  system(paste0("plink --bed ",bed.loc,"ukb_c",i,".bed --bim ",bed.loc,"ukb_c",i,
                ".bim --fam train_split.fam --snps-only just-acgt --make-bed --keep train.10k.fam --extract snps.c",i,
                " --out c",i,".10k"))
  system(paste0("plink --bed ",bed.loc,"ukb_c",i,".bed --bim ",bed.loc,"ukb_c",i,
                ".bim --fam train_split.fam --snps-only just-acgt --keep train.10k.fam",
                " --make-just-bim --extract snps.c",
                i," --out c",i,".cid"))
  
  #load in new bim
  bim = fread(paste0("/home/user/eric/ukbiobank/t1d/data/c",i,".cid.bim"))
  colnames(bim) <- c("chr","snp","cm","start","alt","ref")
  snps <- bim[,2]
  write.table(snps, paste0("snps.c",i), quote = F, col.names = F, row.names = F)
  data <- merge.data.table(bim, snps, by = "snp")
  data <- data[!duplicated(data$snp),]
  data$end <- data$start
  
  
  bim <- data[,-3]
  data$chr <- paste("chr", data$chr, sep = "")
  data <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE,ignore.strand = FALSE)
  
  merged.annotations <- bim
  
  
  # find overlap between data and annotation, hits labeled as 1 and otherwise 0, then reconnect
  wgRna <- fread("~/eric/gcnn/cid_files/wgRna.txt.gz", header=FALSE)
  colnames(wgRna)[2] <- "chr"
  colnames(wgRna)[10] <- "type"
  colnames(wgRna)[3] <- "start"
  colnames(wgRna)[4] <- "end"
  
  mirna <- wgRna[wgRna$type == "miRNA",]
  mirna <- makeGRangesFromDataFrame(mirna, keep.extra.columns = TRUE, ignore.strand = FALSE)
  ol <- findOverlaps(data, mirna)
  hit <- as.data.frame(data[queryHits(ol),])
  if(dim(hit)[1] > 0){
    hit$mirna <- 1
    hit <- hit[,c("snp", "mirna")]
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$mirna <- 0
    nohit <- nohit[,c("snp", "mirna")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  rm(mirna)
  
  cdbox <- wgRna[wgRna$type == "CDBox",]
  cdbox <- makeGRangesFromDataFrame(cdbox, keep.extra.columns = TRUE, ignore.strand = FALSE)
  ol <- findOverlaps(data, cdbox)
  hit <- as.data.frame(data[queryHits(ol),])
  if(dim(hit)[1] > 0){
    hit$cdbox <- 1
    hit <- hit[,c("snp", "cdbox")]
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$cdbox <- 0
    nohit <- nohit[,c("snp", "cdbox")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  rm(cdbox)
  
  hacabox <- wgRna[wgRna$type == "HAcaBox",]
  hacabox <- makeGRangesFromDataFrame(hacabox, keep.extra.columns = TRUE, ignore.strand = FALSE)
  ol <- findOverlaps(data, hacabox)
  hit <- as.data.frame(data[queryHits(ol),])
  if(dim(hit)[1] > 0){
    hit$hacabox <- 1
    hit <- hit[,c("snp", "hacabox")]
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$hacabox <- 0
    nohit <- nohit[,c("snp", "hacabox")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  rm(hacabox)
  
  scarna <- wgRna[wgRna$type == "scaRna",]
  scarna <- makeGRangesFromDataFrame(scarna, keep.extra.columns = TRUE, ignore.strand = FALSE)
  ol <- findOverlaps(data, scarna)
  hit <- as.data.frame(data[queryHits(ol),])
  if(dim(hit)[1] > 0){
    hit$scarna <- 1
    hit <- hit[,c("snp", "scarna")]
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$scarna <- 0
    nohit <- nohit[,c("snp", "scarna")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  rm(scarna, wgRna)
  
  # later might want to change this from "is there overlap" to "what is the highest score"
  dhs <- fread("~/eric/gcnn/cid_files/wgEncodeAwgDnaseMasterSites.txt.gz")
  colnames(dhs)[2] <- "chr"
  colnames(dhs)[6] <- "score"
  colnames(dhs)[3] <- "start"
  colnames(dhs)[4] <- "end"
  dhs <- makeGRangesFromDataFrame(dhs, keep.extra.columns = TRUE, ignore.strand = FALSE)
  ol <- findOverlaps(data, dhs)
  hit <- as.data.frame(data[queryHits(ol),])
  hit$dhs <- 1
  hit <- hit[,c("snp", "dhs")]
  nohit <- as.data.frame(data[-queryHits(ol),])
  nohit$dhs <- 0
  nohit <- nohit[,c("snp", "dhs")]
  rbind.hit <- rbind(hit, nohit)
  merged.annotations <- merge.data.table(merged.annotations, rbind.hit)
  rm(dhs)
  
  test.search <- query(grs, "cpg")
  cpg <- grs[["AH5086"]]
  cpg <- GRanges(cpg)
  ol <- findOverlaps(data, cpg)
  hit <- as.data.frame(data[queryHits(ol),])
  if(dim(hit)[1] > 0){
    hit$cpg <- 1
    hit <- hit[,c("snp", "cpg")]
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$cpg <- 0
    nohit <- nohit[,c("snp", "cpg")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  rm(cpg)
  
  
  
  test.search <- query(grs, "knowngene")
  annot <- grs[["AH5036"]]
  annot <- GRanges(annot)
  ol <- findOverlaps(data, annot)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$knownGene <- 1
    hit <- unique(hit[,c("snp", "knownGene")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$knownGene <- 0
    nohit <- nohit[,c("snp", "knownGene")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, CodingVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$codingVariant <- 1
    hit <- unique(hit[,c("snp", "codingVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$codingVariant <- 0
    nohit <- nohit[,c("snp", "codingVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, by = "snp")
  }
  
  
  
  loc <- locateVariants(data, txdb, IntronVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$IntronVariant <- 1
    hit <- unique(hit[,c("snp", "IntronVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$IntronVariant <- 0
    nohit <- nohit[,c("snp", "IntronVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, FiveUTRVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$FiveUTRVariant <- 1
    hit <- unique(hit[,c("snp", "FiveUTRVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$FiveUTRVariant <- 0
    nohit <- nohit[,c("snp", "FiveUTRVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, ThreeUTRVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$ThreeUTRVariant <- 1
    hit <- unique(hit[,c("snp", "ThreeUTRVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$ThreeUTRVariant <- 0
    nohit <- nohit[,c("snp", "ThreeUTRVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp",all.x = TRUE)
  }
  
  
  
  loc <- locateVariants(data, txdb, SpliceSiteVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$SpliceSiteVariant <- 1
    hit <- unique(hit[,c("snp", "SpliceSiteVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$SpliceSiteVariant <- 0
    nohit <- nohit[,c("snp", "SpliceSiteVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE)
  }
  
  
  
  
  loc <- locateVariants(data, txdb, PromoterVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$PromoterVariant <- 1
    hit <- unique(hit[,c("snp", "PromoterVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$PromoterVariant <- 0
    nohit <- nohit[,c("snp", "PromoterVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, by = "snp",rbind.hit, all.x = TRUE)
  }
  
  
  
  
  test.search <- query(grs, "TFBS")
  annot <- grs[["AH5090"]]
  annot <- GRanges(annot)
  ol <- findOverlaps(data, annot)
  hit <- as.data.frame(data[queryHits(ol),])
  #make this the name of the annotation
  if(dim(hit)[1] > 0){
    hit$TFBS <- 1
    hit <- unique(hit[,c("snp", "TFBS")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$TFBS <- 0
    nohit <- nohit[,c("snp", "TFBS")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE)
  }
  
  merged.annotations <- merged.annotations[!duplicated(merged.annotations),]
  ### end of binary annotations 
  rm(hit, nohit, ol, annot, test.search, rbind.hit, loc, ah, grs)
  
  
  
  
  
  #### CADD
  
  ### fread implementation, current best
  full.merge.cadd <- data.table()
  
  file.list <- list.files(path="/mnt/backup/user/cadd/", pattern = "split.*")
  
  for( k in 1:length(file.list)) {
    print(paste0("Working on section ", k," of ",length(file.list)))
    cadd <- fread(paste0("/mnt/backup/user/cadd/",file.list[k]), showProgress = FALSE)
    colnames(cadd) <- c("chr", "start", "ref", "alt", "rawScore", "cadd")
    #cadd <- subset(cadd, !grepl('^\\d+$', cadd$snp))
    #cadd <- transform(cadd[grep("^\\d+$", cadd$chr),, drop = F], chr = as.numeric(as.character(chr)))
    cadd$chr <- as.integer(cadd$chr)
    merge <- merge.data.table(bim, cadd, by = c("start", "ref", "alt", "chr"))
    bindlist <- list(full.merge.cadd, merge)
    full.merge.cadd <- rbindlist(bindlist)
    rm(cadd, merge, bindlist)
    gc()
  }
  
  
  
  
  merged.annotations <- merge.data.table(merged.annotations, full.merge.cadd, by = c("start","ref","alt","chr","end","snp"), all.x = T)
  merged.annotations <- merged.annotations[,!"rawScore"]
  
  write.table(merged.annotations, paste0("annotations.c",i,".txt"), col.names = T, row.names = F)
  
  
  
  
  #rm(list = setdiff(ls(), c("i", "start_time")))
  gc()
}
