#!/bin/bash

BASE=/hpf/projects/hawkinslab/rob/cna_calls
GENOME=/hpf/projects/hawkinslab/rob/genomes/hg19

mkdir -p $BASE/scripts


cat << EOF > $BASE/scripts/access.sh


#PBS -l vmem=30g,mem=30g
#PBS -l nodes=1:ppn=1
#PBS -o $BASE/logs/access.o
#PBS -j oe
#PBS -N access
set -eu -o pipefail

export PATH=/hpf/projects/hawkinslab/rob/software/miniconda2/bin:\$PATH

cnvkit.py access $GENOME/hg19.fa -o $GENOME/access-10kb-mappable.hg19.bed

EOF


chmod u+x $BASE/scripts/access.sh
COMMAND=$(qsub $BASE/scripts/access.sh)
echo $COMMAND

for i in proton illumina; do
    cat << EOF > $BASE/scripts/$i.cnvkit.sh
#!/bin/bash
#PBS -l vmem=120g,mem=120g
#PBS -l nodes=1:ppn=4
#PBS -l walltime=72:00:00
#PBS -o $BASE/logs/$i.cnvkit.o
#PBS -j oe
#PBS -N $i
# #PBS -W depend=afterok:$COMMAND

set -eu -o pipefail

export PATH=/hpf/projects/hawkinslab/rob/software/miniconda2/bin:\$PATH

cd $BASE/$i/results

cnvkit.py batch $BASE/$i/bam/*tumor.bam \
    --normal $BASE/$i/bam/*normal.bam \
    --targets $BASE/$i/baits.bed \
    --annotate $BASE/refFlat.txt \
    --fasta $GENOME/hg19.fa \
    --drop-low-coverage \
    --access $GENOME/access-10kb-mappable.hg19.bed \
    --output-reference $BASE/$i/results/GRCm38_normal_reference.cnn \
    --output-dir $BASE/$i/results \
    --diagram \
    --scatter \
    -p 4

printf "drawing heatmap\n"
cnvkit.py heatmap -d -o $BASE/$i/heatmap.png $BASE/$i/results/*.cns

printf "calling genders\n"
cd $BASE/$i/results
cnvkit.py sex -o sex.txt *.cnn *.cnr *.cns

printf "done\n"
  
EOF

    chmod u+x $BASE/scripts/$i.cnvkit.sh
    COMMAND2=$(qsub $BASE/scripts/$i.cnvkit.sh)
    echo $COMMAND2

    for SAMPLE in $(cat $BASE/$i/tumors); do
        cat << EOF > $BASE/scripts/$i.$SAMPLE.call.sh
#!/bin/bash

#PBS -l vmem=12g,mem=12g
#PBS -l nodes=1:ppn=1
#PBS -o $BASE/logs/$SAMPLE.$i.calls.o
#PBS -j oe
#PBS -N copycall.$SAMPLE
#PBS -W depend=afterok:33163629

set -eu -o pipefail

export PATH=/hpf/projects/hawkinslab/rob/software/miniconda2/bin:\$PATH

cd $BASE/$i/results

cnvkit.py call -y -m threshold --filter cn -o $SAMPLE.call.cns $SAMPLE.tumor.cns
cnvkit.py gainloss -s $SAMPLE.tumor.cns -o $SAMPLE.gainloss.txt $SAMPLE.tumor.cnr

cnvkit.py scatter -s $SAMPLE.call.cns -o $SAMPLE.scatter.png

echo $SAMPLE copy-calls >> $BASE/$i/cnvkit2.progress

EOF
    

        chmod u+x $BASE/scripts/$i.$SAMPLE.call.sh
        COMMAND3=$(qsub $BASE/scripts/$i.$SAMPLE.call.sh)
        echo $COMMAND3
    done
done



































