library(matrixStats)
library(ggplot2)
library(reshape2)


psiCDF <- function(x, y, z, ...){
  # x is name, y is table, z is type (novel, known, all)
  ppi <- 300
  
  outFile <- paste(x, z, "sdev.cdf.png", sep = ".")
  
  #create input for ggplot2
  
  df.sdev <- melt(data.frame(pHGG = y$T.sdev, Normal = y$N.sdev))
  colnames(df.sdev) <- c("Tissue", "sdev")
  
  sdevTest <- round(ks.test(y$T.sdev, y$N.sdev)$p.value,4)

  plotTitle <- paste(z, x, "eCDF \np =", sdevTest,  sep = " ")
  
  ggplot(df.sdev, aes(x=sdev, color=Tissue)) +
    stat_ecdf(na.rm = T, pad  = F, size = 1) +
    theme_classic() +
    labs(x = "Standard deviation", y = "Density", title = plotTitle) +
    theme(plot.title=element_text( hjust=0.5, vjust=0.5, size = 10, colour="black"),
        axis.text=element_text(size=10, colour="black"),
        axis.title=element_text(size=10, colour="black"),
        aspect.ratio = 1) +
    scale_color_manual(values=c("red", "blue"))
  
  ggsave(outFile, width = 4, height = 4, units = "in", dpi = 300)

}

splicingAmount <- function(x, ...){
   
  # read MATS JCEC output
  matsFile <- paste(x, "MATS.JCEC.txt", sep = ".")
  df.mats <- read.table(matsFile, header = T, row.names = 1)

  # read psi grids
  tumorFile <- paste(x, "tumor.psi", sep = ".")
  df.tumor.psi <- read.table(tumorFile, header = F, sep = '\t')
  normalFile <- paste(x, "normal.psi", sep = ".")
  df.normal.psi <- read.table(normalFile, header = F, sep = '\t')
  
  # create stats tables
  df.stats <- data.frame(ID = df.mats$ID.1,
                         T.mean = rowMeans(df.tumor.psi, na.rm = T),
                         T.median = rowMedians(as.matrix(df.tumor.psi), na.rm = T),
                         T.sdev = rowSds(as.matrix(df.tumor.psi), na.rm = T),
                         N.mean = rowMeans(df.normal.psi, na.rm = T),
                         N.median = rowMedians(as.matrix(df.normal.psi), na.rm = T),
                         N.sdev = rowSds(as.matrix(df.normal.psi), na.rm = T))
  
  diff.mats <- subset(df.mats, FDR < 0.05 & abs(IncLevelDifference) > 0.15)
 
  # draw ecdf curves
  psiCDF(x, df.stats, "all")
   
  return(df.stats)
  
}


SE <- splicingAmount("SE")
MXE <- splicingAmount("MXE")
A5SS <- splicingAmount("A5SS")
A3SS <- splicingAmount("A3SS")
RI <- splicingAmount("RI")


psiCDF("ASE", rbind(SE, MXE, A5SS, A3SS, RI), "all")
psiCDF("ASE", rbind(SE.novel, MXE.novel, RI.novel), "novel")
psiCDF("ASE", rbind(SE.known, MXE.known, A5SS, A3SS, RI.known), "known")               
               
               