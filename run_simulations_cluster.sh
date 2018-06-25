#!/bin/bash

mkdir slurm 2> /dev/null
snakemake -j 2000 \
    --cluster-config cluster.json \
    --local-cores 8 \
    --cluster "{cluster.sbatch} --parsable -p {cluster.partition} --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}" \
    --latency-wait 60 \
    --timestamp \
    --keep-going \
    --restart-times 1 \
    --cluster-status ./cluster_status.py
