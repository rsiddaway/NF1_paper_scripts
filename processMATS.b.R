args <- commandArgs(trailingOnly=TRUE)
options("scipen" = 100)

eventType <- args[1]

load(paste("tmp/", eventType, ".tmp.RData", sep = ""))

if(eventType == "SE"){
  j5a <- read.table(paste("tmp/", eventType, ".5a.scores", sep = ""))
  df$j5a <- j5a$V2

  j3a <- read.table(paste("tmp/", eventType, ".3a.scores", sep = ""))
  df$j3a <- j3a$V2
  
  j5b <- read.table(paste("tmp/", eventType, ".5b.scores", sep = ""))
  df$j5b <- j5b$V2
  
  j3b <- read.table(paste("tmp/", eventType, ".3b.scores", sep = ""))
  df$j3b <- j3b$V2
}

if(eventType == "A3SS"){
  j5a <- read.table(paste("tmp/", eventType, ".5a.scores", sep = ""))
  df$j5a <- j5a$V2
  
  j3a <- read.table(paste("tmp/", eventType, ".3a.scores", sep = ""))
  df$j3a <- j3a$V2
  
  j3b <- read.table(paste("tmp/", eventType, ".3b.scores", sep = ""))
  df$j3b <- j3b$V2
}

if(eventType == "A5SS"){
  j5a <- read.table(paste("tmp/", eventType, ".5a.scores", sep = ""))
  df$j5a <- j5a$V2

  j5b <- read.table(paste("tmp/", eventType, ".5b.scores", sep = ""))
  df$j5b <- j5b$V2
  
  j3a <- read.table(paste("tmp/", eventType, ".3a.scores", sep = ""))
  df$j3a <- j3a$V2
  
}

if(eventType == "MXE"){
  j5a <- read.table(paste("tmp/", eventType, ".5a.scores", sep = ""))
  df$j5a <- j5a$V2
  
  j3a <- read.table(paste("tmp/", eventType, ".3a.scores", sep = ""))
  df$j3a <- j3a$V2
  
  j5b <- read.table(paste("tmp/", eventType, ".5b.scores", sep = ""))
  df$j5b <- j5b$V2
  
  j3b <- read.table(paste("tmp/", eventType, ".3b.scores", sep = ""))
  df$j3b <- j3b$V2
  
  j5c <- read.table(paste("tmp/", eventType, ".5c.scores", sep = ""))
  df$j5c <- j5c$V2
  
  j3c <- read.table(paste("tmp/", eventType, ".3c.scores", sep = ""))
  df$j3c <- j3c$V2
}

if(eventType == "RI"){
  j5a <- read.table(paste("tmp/", eventType, ".5a.scores", sep = ""))
  df$j5a <- j5a$V2
  
  j3a <- read.table(paste("tmp/", eventType, ".3a.scores", sep = ""))
  df$j3a <- j3a$V2
}


#write completed psi table for full and differential results
outputName <- paste(args[1], ".annotated.txt", sep = "")
write.table(df, outputName, quote = FALSE, row.names = FALSE, sep = '\t')
write.table(df.diff, paste('diff', outputName, sep = "."), quote = FALSE, row.names = FALSE, sep = '\t')

save.image(paste("tmp/", eventType, ".tmp.RData", sep = ""))

quit()

