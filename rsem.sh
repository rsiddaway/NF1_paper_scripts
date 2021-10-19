#!/bin/bash
# #PBS -l vmem=64g,mem=64g
# #PBS -l nodes=1:ppn=14
# #PBS -l walltime=24:00:00
#PBS -o $BASE/logs/$i.rsem.o
#PBS -j oe
#PBS -N $i

set -eu -o pipefail

module load picard-tools/2.5.0

PICARD=/hpf/tools/centos6/picard-tools/2.5.0/picard.jar
GENOME=/hpf/projects/hawkinslab/rob/genomes/GRCh37
BASE=/hpf/projects/hawkinslab/rob/rnaseq

cd $BASE/bamOrig/$i

echo fastq
java -Xmx14g \
    -jar $PICARD SamToFastq \
    I=$BASE/bamOrig/$i/$i.bam \
    FASTQ=$BASE/tmp/${i}_R1.fq.gz \
    SECOND_END_FASTQ=$BASE/tmp/${i}_R2.fq.gz \
    UNPAIRED_FASTQ=$BASE/bamOrig/$i/$i.unpaired.fq.gz

echo trimming

java -XX:ParallelGCThreads=1 -Xmx2G -jar /hpf/tools/centos6/mugqic-pipelines/source/resource/software/trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE \
  -threads 12 \
  -phred33 \
  $BASE/tmp/${i}_R1.fq.gz \
  $BASE/tmp/${i}_R2.fq.gz \
  $BASE/raw/$i/$i.trim.pair1.fastq.gz \
  $BASE/raw/$i/$i.trim.single1.fastq.gz \
  $BASE/raw/$i/$i.trim.pair2.fastq.gz \
  $BASE/raw/$i/$i.trim.single2.fastq.gz \
  ILLUMINACLIP:/hpf/tools/centos6/mugqic-pipelines/source/resource/software/trimmomatic/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:15 \
  TRAILING:30 \
  CROP:100 MINLEN:32     
    
rm $BASE/tmp/$i/*.fq.gz

module load rsem
module load samtools
module load bowtie/1.0.1

mkdir -p $BASE/rsem/$i
cd $BASE/rsem/$i

rsem-calculate-expression \
    -p 12 \
    --no-bam-output \
    --estimate-rspd \
    --calc-ci \
    --paired-end \
    --forward-prob 0 \
    --temporary-folder $TMP/$i \
    <(zcat $BASE/raw/$i/$i.trim.pair1.fastq.gz) \
    <(zcat $BASE/raw/$i/$i.trim.pair2.fastq.gz) \
    $GENOME/rsem-ref/grch37.rsem \
    $i
