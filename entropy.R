setwd("/hpf/projects/hawkinslab/rob/rnaseq/DGE_rsem")

library(entropy)
library(edgeR)


iso.fpkm <- read.table('rsem.isoforms.FPKM.assembled', header = T, row.names=1, sep = '\t')

# calculate per-gene Shannon entropy
iso.entropy <- do.call(rbind.data.frame,
             lapply(unique(iso.fpkm$Gene),
                    function(i){
                      outp = apply(iso.fpkm[iso.fpkm$Gene %in% i,2:ncol(iso.fpkm)], 2, entropy);
                      return(outp)}))

colnames(iso.entropy) <- colnames(iso.fpkm[,2:ncol(iso.fpkm)])
row.names(iso.entropy) <- unique(iso.fpkm$Gene)

save.image('entropy.RData')


