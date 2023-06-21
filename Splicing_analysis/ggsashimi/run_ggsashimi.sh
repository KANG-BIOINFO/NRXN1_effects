#!/bin/bash

module load singularity/3.7.0

singularity run -W $PWD -B $PWD:$PWD \
    ggsashimi_latest.sif -b input_bam_eng.tsv \
    -c chr2:50052271-50092494 \
    -O 3 \
    -C 3 \
    -M 5 \
    -A median \
    --width 10 \
    --shrink \
    --fix-y-scale \
    -g gencode.v41.annotation.gtf \
    --ann-height=4 \
    --height=2