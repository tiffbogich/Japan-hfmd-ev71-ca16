#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/ksimplex/$i.json model/theta_kmcmc_cov.json
    cp design_kmcmc_cov.json model/

    cd model
    plom bootstrap theta_kmcmc_cov.json design_kmcmc_cov.json --engine PBS --walltime 72:00:00
    qsub results/kmcmc_cov/design.sh

    cd ../..

done
