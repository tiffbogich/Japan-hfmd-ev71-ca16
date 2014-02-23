#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/kmcmc_cov/$i.json model/theta_kmcmc.json
    cp design_kmcmc.json model/

    cd model
    plom bootstrap theta_kmcmc.json design_kmcmc.json --engine PBS --walltime 72:00:00
    qsub results/kmcmc/design.sh

    cd ../..

done
