args <- commandArgs(trailingOnly=TRUE)
options("scipen" = 100)
eventType <- args[1]
namesFile <- args[2]

# read MATS output file and annotation
matsFileName <- paste(eventType, ".MATS.JCEC.txt", sep = "")
df <- read.table(matsFileName, header = T)
anno <- read.csv(namesFile, header=T, row.names=1)

depth <- read.table("../../depth.txt", header = T, row.names = 1)
depth$sF <- depth$all/50000000

# read other info
novelEventsFileName <- paste("tmp/", eventType, ".novel", sep = "")
IncFileName <- paste("tmp/", eventType, ".Inc", sep = "")
SkipFileName <- paste("tmp/", eventType, ".Skip", sep = "")
PsiFileName <- paste("tmp/", eventType, ".Psi", sep = "")

novelEvents <- read.table(novelEventsFileName, header = F)
Inc <- read.table(IncFileName, header = F, row.names = 1)
Skip <- read.table(SkipFileName, header = F, row.names = 1)
Psi <- read.table(PsiFileName, header = F, row.names = 1)

#add isNovel column
df$isNovel <- ifelse(df$ID %in% novelEvents$V1, 'no', 'yes')
# df$isNovel <- 'no'

#add coverage columns
Coverage <- Inc + Skip
colnames(Coverage) <- paste(anno$sampleID, ".cov", sep = "")
Coverage.norm <- sweep(Coverage, 2, depth$sF, "/")
df <- cbind(df, Coverage.norm)
df$meanCov.T <- rowMeans(Coverage.norm[,1:64], na.rm = TRUE)
df$meanCov.N <- rowMeans(Coverage.norm[,65:84], na.rm = TRUE)
df$meanCov.all <- rowMeans(Coverage.norm, na.rm = TRUE)

#add psi columns
colnames(Psi) <- paste(anno$sampleID, ".psi", sep = "")
df <- cbind(df, Psi)
df$meanPsi.T <- rowMeans(Psi[,1:64], na.rm = TRUE)
df$meanPsi.N <- rowMeans(Psi[,65:84], na.rm = TRUE)
df$meanPsi.all <- rowMeans(Psi, na.rm = TRUE)


#subset for differentially spliced events and add to main table
df.diff <- subset(df, FDR < 0.05, abs(IncLevelDifference) > 0.15)
df$isDiff <- ifelse(df$ID %in% df.diff$ID, 'yes', 'no')
df.diff$isDiff <- 'yes'

save.image(paste("tmp/", eventType, ".tmp.RData", sep = ""))

#make bed files of various junction sets
if(eventType == "SE"){
  S5A <- data.frame(df$chr,
                    df$upstreamEE - 3,
                    df$upstreamEE + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5A, 
              paste("tmp/", eventType, ".5a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S5B <- data.frame(df$chr,
                    df$exonEnd - 3,
                    df$exonEnd + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5B, 
              paste("tmp/", eventType, ".5b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3A <- data.frame(df$chr,
                    df$exonStart_0base - 20,
                    df$exonStart_0base + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3A, 
              paste("tmp/", eventType, ".3a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3B <- data.frame(df$chr,
                    df$downstreamES - 20,
                    df$downstreamES + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3B, 
              paste("tmp/", eventType, ".3b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

if(eventType == "A3SS"){
  S5A <- data.frame(df$chr,
                    df$flankingEE - 3,
                    df$flankingEE + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5A, 
              paste("tmp/", eventType, ".5a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3A <- data.frame(df$chr,
                    df$longExonStart_0base - 20,
                    df$longExonStart_0base + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3A, 
              paste("tmp/", eventType, ".3a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3B <- data.frame(df$chr,
                    df$shortES - 20,
                    df$shortES + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3B, 
              paste("tmp/", eventType, ".3b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

if(eventType == "A5SS"){
  S5A <- data.frame(df$chr,
                    df$longExonEnd - 3,
                    df$longExonEnd + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5A, 
              paste("tmp/", eventType, ".5a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S5B <- data.frame(df$chr,
                    df$shortEE - 3,
                    df$shortEE + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5B, 
              paste("tmp/", eventType, ".5b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3A <- data.frame(df$chr,
                    df$flankingES - 20,
                    df$flankingES + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3A, 
              paste("tmp/", eventType, ".3a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  
}

if(eventType == "MXE"){
  S5A <- data.frame(df$chr,
                    df$upstreamEE - 3,
                    df$upstreamEE + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5A, 
              paste("tmp/", eventType, ".5a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S5B <- data.frame(df$chr,
                    df$Exon1End - 3,
                    df$Exon1End + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5B, 
              paste("tmp/", eventType, ".5b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S5C <- data.frame(df$chr,
                    df$Exon2End - 3,
                    df$Exon2End + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5C, 
              paste("tmp/", eventType, ".5c.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3A <- data.frame(df$chr,
                    df$Exon1Start_0base - 20,
                    df$Exon1Start_0base + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3A, 
              paste("tmp/", eventType, ".3a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3B <- data.frame(df$chr,
                    df$Exon2Start_0base - 20,
                    df$Exon2Start_0base + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3B, 
              paste("tmp/", eventType, ".3b.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  S3C <- data.frame(df$chr,
                    df$downstreamES - 20,
                    df$downstreamES + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3C, 
              paste("tmp/", eventType, ".3c.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}

if(eventType == "RI"){
  S5A <- data.frame(df$chr,
                    df$upstreamEE - 3,
                    df$upstreamEE + 6,
                    df$ID,
                    1,
                    df$strand)
  write.table(S5A, 
              paste("tmp/", eventType, ".5a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  S3A <- data.frame(df$chr,
                    df$downstreamES - 20,
                    df$downstreamES + 3,
                    df$ID,
                    1,
                    df$strand)
  write.table(S3A, 
              paste("tmp/", eventType, ".3a.bed", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

quit()

