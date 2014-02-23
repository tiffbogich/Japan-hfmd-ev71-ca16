#!/bin/bash

cd ../

mkdir -p results/kmcmc_cov

for i in hfmd_*
do

    echo $i
    cd $i/model

    plom reduce design_kmcmc_cov.json -o mle.json
    cp mle.json ../../results/kmcmc_cov/${i}.json
    cd ../../

done
