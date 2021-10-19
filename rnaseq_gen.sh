#!/bin/bash

module load python/2.7.10
module load mugqic-pipelines/latest_hpf

# Total steps: 1-24
/hpf/tools/centos6/mugqic-pipelines/latest_hpf/pipelines/rnaseq/rnaseq.py \
 -s 1-14,22 \
 -o $OUTPUTDIR \
 -j pbs \
 -l debug \
 -d sample.design \
 -r sample.readset \
 -c rnaseq.ini 1> qsub_132.sh 2> 132.log

chmod +x qsub_132.sh
