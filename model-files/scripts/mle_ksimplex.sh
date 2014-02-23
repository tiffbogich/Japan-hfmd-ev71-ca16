#!/bin/bash

cd ../

mkdir -p results/ksimplex

for i in hfmd_*
do

    echo $i
    cd $i/model
    tail -1 results/ksimplex/trace_0.csv | awk 'BEGIN { FS=","; OFS="\n" } { print $NF }'

    plom pipe results/ksimplex/theta.json -T results/ksimplex/trace_0.csv -o mle.json
    cp mle.json ../../results/ksimplex/${i}.json
    cd ../../

done
