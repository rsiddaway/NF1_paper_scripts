#!/bin/bash
#PBS -o /hpf/projects/hawkinslab/rob/rnaseq/logs/assemble.o
#PBS -j oe
#PBS -l vmem=10g,mem=10g
#PBS -N rsem.assemble

set -eu -o pipefail

BASE=/hpf/projects/hawkinslab/rob/rnaseq

cd /hpf/projects/hawkinslab/rob/rnaseq/DGE_rsem
mkdir -p tmp


for j in isoforms; do
    # for j in isoforms genes; do
    for k in TPM FPKM; do
        cat rsem.$j.IDs > rsem.$j.$k.assembled
    done
done

for i in $(cat $BASE/final_samples); do
    for j in genes isoforms; do
        for k in TPM FPKM; do
            
            if [ "$k" == "TPM" ]; then
                cutCol=-f6
            else cutCol=-f7
            fi
            
            printf "$i\n" > tmp/$i.$j.$k
            
            tail -n+2 /hpf/projects/hawkinslab/rob/rnaseq/rsem/$i/$i.$j.results | cut $cutCol >> tmp/$i.$j.$k
            paste rsem.$j.$k.assembled tmp/$i.$j.$k > tmp/$j.$k
            cat tmp/$j.$k > rsem.$j.$k.assembled
            
        done
        
        if [ "$j" == "isoforms" ]; then
            k=IsoPct
            cat rsem.$j.IDs > rsem.$j.$k.assembled
            cutCol=-f8
            printf "$i\n" > tmp/$i.$j.$k
            
            tail -n+2 /hpf/projects/hawkinslab/rob/rnaseq/rsem/$i/$i.$j.results | cut $cutCol >> tmp/$i.$j.$k
            paste rsem.$j.$k.assembled tmp/$i.$j.$k > tmp/$j.$k
            cat tmp/$j.$k > rsem.$j.$k.assembled
        fi
        
    done
done

cat rsem.isoforms.IDs > rsem.isoforms.IsoPct.assembled
k=IsoPct
for i in $(cat $BASE/final_samples); do
    for j in isoforms; do
                
        printf "$i\n" > tmp/$i.$j.$k
        
        tail -n+2 /hpf/projects/hawkinslab/rob/rnaseq/rsem/$i/$i.$j.results | cut -f8 >> tmp/$i.$j.IsoPct
        paste rsem.$j.$k.assembled tmp/$i.$j.$k > tmp/$j.$k
        
        cat tmp/$j.$k > rsem.$j.$k.assembled
        
    done
done


rm -rf tmp

# module load R
# Rscript entropy.R


