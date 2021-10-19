#!/bin/bash

cd /hpf/projects/hawkinslab/rob/rnaseq/rmats/hgg_vs_n
chmod u+x diff.cosmic.r

for i in cat SE A5SS A3SS MXE RI; do

  cat << EOF > tmp/$i.cosmic.sh
#!/bin/bash
#PBS -l vmem=4g,mem=4g
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=1
#PBS -N $i
#PBS -o /hpf/projects/hawkinslab/rob/rnaseq/rmats/hgg_vs_n/tmp/$i.cosmic.o
#PBS -j oe

cd /hpf/projects/hawkinslab/rob/rnaseq/rmats/hgg_vs_n
module load R
Rscript diff.cosmic.r $i

EOF

    chmod u+x tmp/$i.cosmic.sh
    qsub tmp/$i.cosmic.sh

done

