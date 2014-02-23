#!/bin/bash

cd ../

mkdir -p results/pmcmc_ode_cov

for i in hfmd_*
do

    echo $i
    cd $i/model

    plom reduce design_pmcmc_ode_cov.json -o mle.json
    cp mle.json ../../results/pmcmc_ode_cov/${i}.json
    cd ../../

done
