ppi = 300
library(pheatmap)
library(stringr)
library(maser)
library(rtracklayer)

hgg <- maser(getwd(), c("pHGG", "Normal"), ftype = "JCEC")

ens_gtfb <- rtracklayer::import.gff('Homo_sapiens.GRCh38.99.gtf')

x.filt<-geneEvents()

plot.maser <- function(x, v, y, z, ...){
  #v: gene
  #x: maser object
  #y: event ID
  #z: event cat i.e. SE etc
  x.filt <- geneEvents(v, x)
  x.mapped <- mapTranscriptsToEvents(x.filt, ens_gtfb)
  plotUniprotKBFeatures(x.mapped,
  z,
  event_id = y,
  gtf = ens_gtfb,
  features = c("mod-res"),
  show_transcripts = TRUE)
  # features = c("domain", "region"),
}


#run separately for mod-res and domain/region
pdf('chromatin.pdf')
    #insert events.full.maser
dev.off()
