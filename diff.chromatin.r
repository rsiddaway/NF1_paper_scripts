args <- commandArgs(trailingOnly=TRUE)

eventType <- args[1]

chromatinGenes <- data.frame(read.delim('chromatin_genes', header=F))

inFile <- paste("diff", eventType, "annotated.txt", sep = ".")

df <- read.delim(inFile, header = T, sep='\t')

df.chromatin <- df[df$geneSymbol %in% chromatinGenes$V1,]

df.chromatin.up <- subset(df.chromatin, IncLevelDifference > 0.15)
nrow(df.chromatin.up)

df.chromatin.down <- subset(df.chromatin, IncLevelDifference < -0.15)
nrow(df.chromatin.down)

df.chromatin.out <- rbind(df.chromatin.up, df.chromatin.down)
nrow(df.chromatin.out)

outputName <- paste("chromatin", eventType, "annotated.txt", sep = ".")
write.table(df.chromatin.out, outputName, quote = FALSE, row.names = FALSE, sep = '\t')

# write full table for processing with maser
inFile.b <- paste(eventType, ".MATS.JCEC.txt", sep = "")
df.b <- data.frame(read.delim(inFile.b, header = T, sep='\t'))
df.b.chromatin <- df.b[df.b$ID %in% df.chromatin.out$ID,]

outputName.b <- paste("chromatin", eventType, "MATS.JCEC.txt", sep = ".")
# write.table(df.b.chromatin, outputName.b, quote = FALSE, row.names = FALSE, sep = '\t')
write.table(df.b.chromatin, outputName.b, row.names = FALSE, sep = '\t')

quit()

