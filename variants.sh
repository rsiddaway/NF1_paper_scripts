#!/bin/bash

BASE=$(pwd)

mkdir -p $BASE/tmp
TMPDIR=$BASE/tmp

mkdir -p $BASE/hc
HC=$BASE/hc

GENOME=/path/to/Homo_sapiens.hg19.fa
DBSNP=/hpf/tools/centos6/mugqic-pipelines/source/resource/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.dbSNP144.vcf
SNPEFF_HOME=/hpf/tools/centos6/mugqic-pipelines/source/resource/software/snpEff/snpEff_4_2


gatk HaplotypeCaller \
    --sample-name $i \
    -R $GENOME/Homo_sapiens.hg19.fa \
    -D $DBSNP \
    -I /localhd/\$PBS_JOBID/$i.bam \
    -O $i.hc.vcf.gz \
    --TMP_DIR /localhd/\$PBS_JOBID \
    --read-validation-stringency SILENT
    
touch filtering.$i

gatk SelectVariants \
    -R $GENOME/Homo_sapiens.hg19.fa \
    -V $i.hc.vcf.gz \
    --select-type-to-include SNP \
    -O $i.snp.vcf.gz

gatk SelectVariants \
    -R $GENOME/Homo_sapiens.hg19.fa \
    -V $i.hc.vcf.gz \
    --select-type-to-include INDEL \
    -O $i.indel.vcf.gz

gatk SelectVariants \
    -R $GENOME/Homo_sapiens.hg19.fa \
    -V $i.hc.vcf.gz \
    --select-type-to-include MIXED \
    -O $i.mixed.vcf.gz

gatk VariantFiltration \
    -V $i.snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -filter "DP > 20" --filter-name "depth-20" \
    -O $i.snps_filtered.vcf.gz

gatk VariantFiltration \
    -V $i.indel.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -filter "DP > 20" --filter-name "depth-20" \
    -O $i.indels_filtered.vcf.gz

gatk VariantFiltration \
    -V $i.mixed.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -filter "DP > 20" --filter-name "depth-20" \
    -O $i.mixed_filtered.vcf.gz


java -jar /hpf/tools/centos6/picard-tools/2.5.0/picard.jar MergeVcfs \
    I=$i.snps_filtered.vcf.gz \
    I=$i.indels_filtered.vcf.gz \
    I=$i.mixed_filtered.vcf.gz \
    O=$i.variants.vcf.gz \
	SEQUENCE_DICTIONARY=$GENOME/Homo_sapiens.hg19.dict


java -jar -Xmx10G /hpf/tools/centos6/varscan/2.3.6/VarScan.jar somatic \
    $indir/$normal $indir/$tumor $outdir/$out

####
samtools mpileup -B -f $GENOME /localhd/\$PBS_JOBID/$i.bam | \
    java -Xmx10G /hpf/tools/centos6/varscan/2.3.6/VarScan.jar mpileup2snp --output-vcf 1

###
java -Djava.io.tmpdir=/localhd/\$PBS_JOBID -XX:ParallelGCThreads=2 -Xmx8G \
  -jar $SNPEFF_HOME/SnpSift.jar annotate \
  $DBSNP \
  $i.$caller.variants.vcf \
  > $i.$caller.variants.snpeff.vcf
