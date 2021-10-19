#!/bin/bash
#PBS -l vmem=100g,mem=100g
#PBS -l nodes=1:ppn=16
#PBS -o $BASE/rmats/$contr/rmats.o
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -N rmats

set -eu -o pipefail

module load rMATS
module load R/3.5.1
module load bedtools

BASE=/hpf/projects/hawkinslab/rob/rnaseq
GENOME=/hpf/projects/hawkinslab/rob/genomes/GRCh37
GENOMEFASTA=$GENOME/Homo_sapiens.GRCh37.fa

mkdir -p $BASE/rmats

cd $BASE/rmats

#$sampA and $sampB are comma-separated lists of bam files

python /hpf/tools/centos6/rMATS/4.0.1/rmats.py \
    --b1 $sampA \
    --b2 $sampB \
    --gtf $GENOME/Homo_sapiens.GRCh37.75.gtf \
    --od $BASE/rmats \
    -t paired \
    --libType fr-firststrand \
    --readLength 101 \
    --nthread 16 \
    --tstat 6

module purge
module load bedtools
module load R

mkdir -p $BASE/rmats/tmp

for event in SE MXE A3SS A5SS RI; do

    echo $event
    
    #cut novel events
    tail -n+2 fromGTF.novelEvents.$event.txt | cut -f1,2 > tmp/$event.novel
    
    #cut coverages and psi for tumors and normals
    if [ "$event" == "MXE" ]; then
        Inc=15,17
        Skip=16,18
        Psi=23,24
    else
        Inc=13,15
        Skip=14,16
        Psi=21,22
    fi
    
    tail -n+2 $event.MATS.JCEC.txt | cut -f1,$Inc | tr "," "\t" > tmp/$event.Inc
    tail -n+2 $event.MATS.JCEC.txt | cut -f1,$Skip | tr "," "\t" > tmp/$event.Skip
    tail -n+2 $event.MATS.JCEC.txt | cut -f1,$Psi | tr "," "\t" > tmp/$event.Psi

    
    # run processing steps
    echo processing
    Rscript processMATS.R $event anno.csv
    
    # calculate splice site strengths with maxEntScan
    ls tmp/$event.5*.bed | cut -f2 -d '/' | sed 's/.bed//' > tmp/5ss.$event.txt
    ls tmp/$event.3*.bed | cut -f2 -d '/' | sed 's/.bed//' > tmp/3ss.$event.txt
    echo junctions
    for juncs in $(cat tmp/5ss.$event.txt); do
        sed 's/chr//' tmp/$juncs.bed | \
        bedtools getfasta -fi $GENOMEFASTA -bed - -s -tab | \
        cut -f2 > tmp/$juncs.seqs

        echo scores
        cd /hpf/projects/hawkinslab/rob/software/maxentscan
        perl score5.pl $BASE/tmp/$juncs.seqs > $BASE/tmp/$juncs.scores
        
        cd $BASE
    done
    for juncs in $(cat tmp/3ss.$event.txt); do
        echo fasta
        sed 's/chr//' tmp/$juncs.bed | \
        bedtools getfasta -fi $GENOMEFASTA -bed - -s -tab | \
        cut -f2 > tmp/$juncs.seqs

        echo scores
        cd /hpf/projects/hawkinslab/rob/software/maxentscan
        perl score3.pl $BASE/tmp/$juncs.seqs > $BASE/tmp/$juncs.scores
        
        cd $BASE
    done
    echo final annotation
    # add splice sites to annotated data table
    Rscript processMATS.b.R $event

# done

echo combining event types

Rscript consolidate.R SE.annotated.txt A3SS.annotated.txt A5SS.annotated.txt MXE.annotated.txt RI.annotated.txt anno.csv

echo finished

# to run maser on subset of events
# ./lift.sh

# to test AS stdev
# Rscript psiVariability.R


