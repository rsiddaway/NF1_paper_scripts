args <- commandArgs(trailingOnly=TRUE)

options("scipen" = 100)

SE.file <- args[1]
A3SS.file <- args[2]
A5SS.file <- args[3]
MXE.file <- args[4]
RI.file <- args[5]
namesFile <- args[6]

SE <- read.table(SE.file, header = T)
A3SS <- read.table(A3SS.file, header = T)
A5SS <- read.table(A5SS.file, header = T)
MXE <- read.table(MXE.file, header = T)
RI <- read.table(RI.file, header = T)

anno <- read.csv(namesFile, header=T, row.names=1)

depth <- read.table("../../depth.txt", header = T, row.names = 1)
depth$sF <- depth$all/50000000

rm(SE.file, A3SS.file, A5SS.file, MXE.file, RI.file, namesFile)

save.image("processed.RData")

quit()
