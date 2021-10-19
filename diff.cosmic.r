args <- commandArgs(trailingOnly=TRUE)

eventType <- args[1]

cosmicGenes <- read.delim('Census_allMon_Apr_13_02_55_32_2020.tsv', header=T, sep='\t')

inFile <- paste("diff", eventType, "annotated.txt", sep = ".")

df <- read.delim(inFile, header = T, sep='\t')

df.cosmic <- df[df$geneSymbol %in% cosmicGenes$GeneSymbol,]

df.cosmic.up <- subset(df.cosmic, IncLevelDifference > 0.15)
nrow(df.cosmic.up)

df.cosmic.down <- subset(df.cosmic, IncLevelDifference < -0.15)
nrow(df.cosmic.down)

df.cosmic.out <- rbind(df.cosmic.up, df.cosmic.down)
nrow(df.cosmic.out)

outputName <- paste("cosmic", eventType, "annotated.txt", sep = ".")
write.table(df.cosmic.out, outputName, quote = FALSE, row.names = FALSE, sep = '\t')

# write full table for processing with maser
inFile.b <- paste(eventType, ".MATS.JCEC.txt", sep = "")
df.b <- data.frame(read.delim(inFile, header = T, sep='\t'))
df.b.cosmic <- df.b[df.b$ID %in% df.cosmic.out$ID,]

outputName.b <- paste("cosmic", eventType, "MATS.JCEC.txt", sep = ".")
write.table(df.b.cosmic, outputName.b, quote = FALSE, row.names = FALSE, sep = '\t')

quit()

