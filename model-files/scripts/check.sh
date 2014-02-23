#!/bin/bash

for i in hfmd_*
do

    echo $i
    echo

    for j in `seq 0 4`;
    do
        echo $i $j
        tail -1 ${i}/model/results/kmcmc_cov/covariance_$j.csv
        echo ---
        tail -1 ${i}/model/results/pmcmc_ode_cov/covariance_$j.csv
        echo ---
    done


done
