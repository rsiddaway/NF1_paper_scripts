#!/bin/bash
#PBS -l vmem=10g,mem=10g
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=1

set -eu -o pipefail

module load R

cd /hpf/projects/hawkinslab/rob/rnaseq/rmats

for i in SE A3SS A5SS RI MXE; do
    echo $i
    Rscript diff.chromatin.r $i
done

#####
#change to cosmic for cosmic genes
#####
mkdir -p chromatin
cd chromatin


for i in SE A3SS A5SS RI; do

    mv ../chromatin.$i.MATS.JCEC.txt tmp.$i
    
    #extract coordinates
    cut -f4,6,7,12 tmp.$i | tail -n+2 > $i.1.bed
    cut -f4,8,9,12 tmp.$i | tail -n+2 > $i.2.bed
    cut -f4,10,11,12 tmp.$i | tail -n+2 > $i.3.bed
    
    for j in 1 2 3; do
        ./liftOver $i.$j.bed hg19ToHg38.over.chain $i.$j.lift.bed $i.$j.unmapped
        cut -f2,3 $i.$j.lift.bed > $i.$j.B
    done
    
    cut -f1-5 tmp.$i > $i.A
    head -n 1 tmp.$i | cut -f6-11 > $i.B
    paste $i.1.B $i.2.B $i.3.B >> $i.B
    cut -f12-23 tmp.$i > $i.C
    
    paste $i.A $i.B $i.C > $i.MATS.JCEC.txt

done

for i in MXE; do

    mv ../chromatin.$i.MATS.JCEC.txt tmp.$i
    
    #extract coordinates
    cut -f4,6,7,14 tmp.$i | tail -n+2 > $i.1.bed
    cut -f4,8,9,14 tmp.$i | tail -n+2 > $i.2.bed
    cut -f4,10,11,14 tmp.$i | tail -n+2 > $i.3.bed
    cut -f4,12,13,14 tmp.$i | tail -n+2 > $i.4.bed
    
    for j in 1 2 3 4; do
        ./liftOver $i.$j.bed hg19ToHg38.over.chain $i.$j.lift.bed $i.$j.unmapped
        cut -f2,3 $i.$j.lift.bed > $i.$j.B
    done
    
    cut -f1-5 tmp.$i > $i.A
    head -n 1 tmp.$i | cut -f6-13 > $i.B
    paste $i.1.B $i.2.B $i.3.B $i.4.B >> $i.B
    cut -f14-25 tmp.$i > $i.C
    
    paste $i.A $i.B $i.C > $i.MATS.JCEC.txt

done
for i in RI SE MXE RI A3SS A5SS; do
    tail -n+2 $i.MATS.JCEC.txt | awk -v r=$i -v OFS='\t' '{print $3, "hgg", $1, r}' >> events.maser
done

IFS='
'
for i in $(cat events.maser); do
gene=$(echo $i | cut -f1)
cate=$(echo $i | cut -f2)
evid=$(echo $i | cut -f3)
evtype=$(echo $i | cut -f4)

printf "plot.maser('$gene',$cate,$evid,'$evtype')\n" >> events.full.maser
done    
    
#add events into runMaser.R for analysis
#Rscript runMaser.R    
