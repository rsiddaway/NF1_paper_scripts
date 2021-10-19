#!/bin/bash

BASE=/hpf/projects/hawkinslab/rob/rnaseq
CHIP=/hpf/projects/hawkinslab/rob/software/gsea
GMT=/hpf/projects/hawkinslab/rob/software/gsea/msigdb_v6.0_GMTs

for i in logs scripts; do
    rm -rf $BASE/$i/gsea
    mkdir -p $BASE/$i/gsea
done

cd $BASE


for i in $(cat contrasts); do
    
    cat << EOF > $BASE/scripts/gsea/$i.rnk.r
#!/usr/bin/env R
setwd("$BASE/DGE/$i")

library(data.table)

df.in <- read.table("dge_results.csv", sep = '\t', header = T)
df <- (setDT(df.in)[, .SD[which.max(log_CPM)], by = gene_symbol])


df.A <- data.frame(cbind(as.character(df\$gene_symbol), -log10(df\$edger.adj.p.value) * sign(df\$log_FC)))
write.table(df.A, '$i.log10.rnk', sep='\t', quote=FALSE, col.names = FALSE, row.names = FALSE)


EOF

    cat << EOF > $BASE/scripts/gsea/$i.rnk.sh
#!/bin/bash
#PBS -l vmem=20g,mem=20g
#PBS -l nodes=1:ppn=1
#PBS -N $i.rnk
#PBS -o $BASE/logs/gsea/$i.rnk.o
#PBS -j oe

module load R/3.4.2

cd $BASE
chmod u+x $BASE/scripts/gsea/$i.rnk.r
Rscript $BASE/scripts/gsea/$i.rnk.r

EOF

    chmod u+x $BASE/scripts/gsea/$i.rnk.sh
    COMMAND=$(qsub $BASE/scripts/gsea/$i.rnk.sh)
    echo $COMMAND

    for k in log10; do
        for j in $(cat pathway_sets); do
            mkdir -p $BASE/DGE/$i/gsea/$k/$j
        
        cat << EOF > $BASE/DGE/$i/gsea/$i.$j.$k.gsea.sh
#!/bin/bash
#PBS -l vmem=64g,mem=64g
#PBS -l nodes=1:ppn=1
#PBS -l walltime=16:00:00
#PBS -N $i.$j.$k
#PBS -o $BASE/DGE/$i/gsea/$i.$j.$k.o
#PBS -j oe
# #PBS -W depend=afterok:$COMMAND

set -eu -o pipefail

module load gsea2/2.2.3

cd $BASE/DGE/$i/gsea/$k/$j

java -Xmx48G \
-cp /hpf/tools/centos6/gsea2/2.2.3/gsea2-2.2.3.jar \
xtools.gsea.GseaPreranked \
-rnk $BASE/DGE/$i/$i.$k.rnk \
-gmx $GMT/$j.v6.0.symbols.gmt \
-chip $CHIP/ENSEMBL_human_gene.chip \
-collapse false \
-norm meandiv \
-scoring_scheme classic \
-make_sets true \
-mode Max_probe \
-gui false \
-rpt_label $i.$j.$k \
-out $BASE/DGE/$i/gsea/$k/$j \
-set_min 15 \
-set_max 500 \
-nperm 1000 \
-rnd_seed 123 \
-zip_report false \
-include_only_symbols true \
-plot_top_x 1000
echo finished
EOF


        chmod u+x $BASE/DGE/$i/gsea/$i.$j.$k.gsea.sh
        qsub $BASE/DGE/$i/gsea/$i.$j.$k.gsea.sh
        done
        
    done
done

