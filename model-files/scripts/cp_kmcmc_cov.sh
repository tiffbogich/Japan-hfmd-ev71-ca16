#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/ksimplex/$i.json model/theta_kmcmc.json
    cp design_kmcmc.json model/

    cd model

    fit bootstrap theta_kmcmc.json design_kmcmc.json --queue PBS,71:00:00
    mv run_designs.sh run_kmcmc.sh

    cp results/kmcmc_cov/*.csv results/kmcmc/

    cd ../..

done
