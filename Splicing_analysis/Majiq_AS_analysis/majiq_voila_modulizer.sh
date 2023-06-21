#!/bin/bash

#BSUB -M 128000
#BSUB -W 12:00
#BSUB -n 12
#BSUB -R "rusage[mem=128000] span[hosts=1]"
#BSUB -J run_voila_modulizer
#BSUB -o logging/run_voila_modulizer.out
#BSUB -e logging/run_voila_modulizer.err

module load anaconda3
source activate /data/aronow/Kang/py_envir/splice/

voila modulize \
    'output/spliceGraph/splicegraph.sql' \
    'output/heterogen/A_HET-A_WT.het.voila' \
    -d 'output/voila_modulized' --overwrite