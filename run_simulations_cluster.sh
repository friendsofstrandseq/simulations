#!/bin/bash

snakemake -j 1000 \
    --cluster-config cluster.json
    --local-cores 1 \
    --cluster "{cluster.sbatch} -p {cluster.partition} --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}" \
    --latency-wait 60 \
    --timestamp \
    --keep-going \
    --restart-times 1
