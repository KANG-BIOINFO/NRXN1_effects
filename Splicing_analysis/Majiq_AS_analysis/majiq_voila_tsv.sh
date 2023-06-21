#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 12
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J run_voila_tsv
#BSUB -o logging/run_voila_tsv.out
#BSUB -e logging/run_voila_tsv.err

module load anaconda3
source activate /data/aronow/Kang/py_envir/splice/

voila tsv \
    output/spliceGraph/splicegraph.sql \
    output/heterogen/*.het.voila \
    --show-all \
    -f output/voila_tsv/all-pairs-het.tsv