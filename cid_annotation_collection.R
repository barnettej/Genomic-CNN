
#library(readr)
library(AnnotationHub)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(data.table)
library(plyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(data.table)
library(doParallel)
library(BBmisc)

system("export PATH=$PATH:/home/user/user/bin")
run_name <- "t2d"
gene_annot_file <- paste0(run_name,".magma.genes.annot")
# gene_loc <- "/home/user/user/gsprs/req_files/NCBI38.gene.loc"
gene_loc <- "/home/user/user/gsprs/req_files/NCBI37.3.gene.loc"
gene_set <- "/home/user/user/gsprs/req_files/c5.all.v7.1.entrez.gmt"
# gene_set <- "/home/user/user/gsprs/req_files/h.all.v2022.1.Hs.entrez.gmt"
magma_bfile <- paste0("/home/user/user/ukbiobank/t2d/cid/",run_name,".train")
ah <- AnnotationHub()
grs <- ah[ah$rdataclass == "GRanges" & ah$genome == "hg19",]
number_of_sets <- 20

# Load hg19-compatible dataset
setwd("/home/user/user/ukbiobank/t2d/data")
system("mkdir plink_files")
bed.loc <- "/mnt/backup/user/gfetch/"

#txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

for(i in 22:1){
  ah <- AnnotationHub()
  grs <- ah[ah$rdataclass == "GRanges" & ah$genome == "hg19",]
  
  
  # Load hg19-compatible dataset
  setwd("/home/user/user/ukbiobank/t2d/data")
  bed.loc <- "/mnt/backup/user/gfetch/"
  
  #txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  
  start_time <- Sys.time()
  data = fread(paste0("/home/user/user/ukbiobank/imputation_info/ukb_mfi_chr",i,"_v3.txt"))
  bim = fread(paste0(bed.loc,"ukb_c",i,".bim"))
  colnames(bim) <- c("chr","snp","cm","start","alt","ref")
  colnames(data) <- c("change","snp","start", "ref","alt","maf","major","impute_score")
  
  

  #load in summary statistics file so we can match with our data
  sum_stat <- fread("/home/user/user/ukbiobank/t2d/cid/t2d.sumstat")
  colnames(sum_stat)[1] <- "snp"

  #merge and keep only snps that are found in both
  # data <- merge(data, sum_stat, by = "snp")
  snps <- data[,2]
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
  bim = fread(paste0("/home/user/user/ukbiobank/t2d/data/c",i,".cid.bim"))
  colnames(bim) <- c("chr","snp","cm","start","alt","ref")
  snps <- bim[,2]
  write.table(snps, paste0("snps.c",i), quote = F, col.names = F, row.names = F)
  data <- merge.data.table(bim, snps, by = "snp", sort = F)
  data <- data[!duplicated(data$snp),]
  data$end <- data$start
  
  
  bim <- data[,-3]
  data$chr <- paste("chr", data$chr, sep = "")
  data <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE,ignore.strand = FALSE)
  
  merged.annotations <- bim
  
  
  # find overlap between data and annotation, hits labeled as 1 and otherwise 0, then reconnect
  wgRna <- fread("~/user/gcnn/cid_files/wgRna.txt.gz", header=FALSE)
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
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort = F)
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
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort =F)
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
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort = F)
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
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort = F)
  }
  rm(scarna, wgRna)
  

  dhs <- fread("~/user/gcnn/cid_files/wgEncodeAwgDnaseMasterSites.txt.gz")
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
  merged.annotations <- merge.data.table(merged.annotations, rbind.hit, sort = F)
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
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort = F)
  }
  rm(cpg)
  
  
  
  test.search <- query(grs, "knowngene")
  annot <- grs[["AH5036"]]
  annot <- GRanges(annot)
  ol <- findOverlaps(data, annot)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$knownGene <- 1
    hit <- unique(hit[,c("snp", "knownGene")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$knownGene <- 0
    nohit <- nohit[,c("snp", "knownGene")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, sort = F)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, CodingVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$codingVariant <- 1
    hit <- unique(hit[,c("snp", "codingVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$codingVariant <- 0
    nohit <- nohit[,c("snp", "codingVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, all.x = TRUE, by = "snp", sort = F)
  }
  
  
  
  loc <- locateVariants(data, txdb, IntronVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$IntronVariant <- 1
    hit <- unique(hit[,c("snp", "IntronVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$IntronVariant <- 0
    nohit <- nohit[,c("snp", "IntronVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE, sort = F)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, FiveUTRVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$FiveUTRVariant <- 1
    hit <- unique(hit[,c("snp", "FiveUTRVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$FiveUTRVariant <- 0
    nohit <- nohit[,c("snp", "FiveUTRVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE, sort = F)
  }
  
  
  
  
  
  loc <- locateVariants(data, txdb, ThreeUTRVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$ThreeUTRVariant <- 1
    hit <- unique(hit[,c("snp", "ThreeUTRVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$ThreeUTRVariant <- 0
    nohit <- nohit[,c("snp", "ThreeUTRVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp",all.x = TRUE, sort = F)
  }
  
  
  
  loc <- locateVariants(data, txdb, SpliceSiteVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$SpliceSiteVariant <- 1
    hit <- unique(hit[,c("snp", "SpliceSiteVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$SpliceSiteVariant <- 0
    nohit <- nohit[,c("snp", "SpliceSiteVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE, sort = F)
  }
  
  
  
  
  loc <- locateVariants(data, txdb, PromoterVariants())
  ol <- findOverlaps(data, loc)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$PromoterVariant <- 1
    hit <- unique(hit[,c("snp", "PromoterVariant")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$PromoterVariant <- 0
    nohit <- nohit[,c("snp", "PromoterVariant")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, by = "snp",rbind.hit, all.x = TRUE, sort = F)
  }
  
  
  
  
  test.search <- query(grs, "TFBS")
  annot <- grs[["AH5090"]]
  annot <- GRanges(annot)
  ol <- findOverlaps(data, annot)
  hit <- as.data.frame(data[queryHits(ol),])

  if(dim(hit)[1] > 0){
    hit$TFBS <- 1
    hit <- unique(hit[,c("snp", "TFBS")])
    nohit <- as.data.frame(data[-queryHits(ol),])
    nohit$TFBS <- 0
    nohit <- nohit[,c("snp", "TFBS")]
    rbind.hit <- rbind(hit, nohit)
    merged.annotations <- merge.data.table(merged.annotations, rbind.hit, by = "snp", all.x = TRUE, sort = F)
  }
  
  merged.annotations <- merged.annotations[!duplicated(merged.annotations),]
  ### end of binary annotations 
  rm(hit, nohit, ol, annot, test.search, rbind.hit, loc, ah, grs)
  
  
  magma_out <- fread(paste0(run_name,".magma.sets.out"),skip = 3)
  magma_out <- magma_out[order(SE, decreasing = F)]
  magma_cut <- magma_out[1:(dim(magma_out)[1]/10),]
  magma_cut <- magma_cut[order(P)]
  magma_cut <- magma_cut[1:number_of_sets,]
  #make list of snps for plink command
  gene_sets_snps <- fread(paste0(run_name,".gene_sets_snps.annot"), header = FALSE)
  colnames(gene_sets_snps) <- c("FULL_NAME", "SNPs")
  merged <- merge(magma_cut, gene_sets_snps, by.x = "FULL_NAME", by.y = "FULL_NAME", sort = F)
  plink_sets <- data.table(merged)
  plink_sets$space_SNPs <- gsub("\\t"," ", plink_sets$SNPs)
  plink_sets$space_SNPs <- gsub("\\s+", " ", plink_sets$space_SNPs)
  plink_sets$unique_SNPs <- vapply(lapply(strsplit(plink_sets$space_SNPs, " "), unique), paste, character(1L), collapse = " ")
  
  ### gene set annotation
  for(j in 1:number_of_sets) {
    set_name = plink_sets$FULL_NAME[j]
    write.table(plink_sets$unique_SNPs[j], paste("plink_files/snps_",set_name, sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
    system(paste("tr \" \" \"\n\" < plink_files/snps_",set_name," > plink_files/set_",set_name,".txt", sep = ""))
    system(paste("rm plink_files/snps_",set_name,sep = ""))
    
    snps_in_set <- fread(paste0("plink_files/set_",set_name,".txt"),header = F)
    set.annot <- data.table()
    set.annot$snp <- bim$snp
    set.annot$annot <- as.integer(set.annot$snp %in% snps_in_set$V1)
    colnames(set.annot)[2] <- set_name
    summary(set.annot)
    merged.annotations <- merge.data.table(merged.annotations, set.annot, by = "snp", all.x = TRUE, sort = F)
  }
  
  write.table(merged.annotations, paste0(run_name,".annotations.c",i,".txt"), col.names = T, row.names = F)
  

  gc()
}
