#!/bin/bash
set -x

BATCH_SIZE=5

nextflow run main.nf \
    --input data/samples_to_run.csv \
    --num_samples $BATCH_SIZE \
    -profile docker \
    -entry fetch_sra &
    
BACK_PID=$!
wait $BACK_PID
rm -rf work


