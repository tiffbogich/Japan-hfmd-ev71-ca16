#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/pmcmc_ode_cov/$i.json model/theta_pmcmc_ode.json
    cp design_pmcmc_ode.json model/

    cd model
    plom bootstrap theta_pmcmc_ode.json design_pmcmc_ode.json --engine PBS --walltime 72:00:00
    qsub results/pmcmc_ode/design.sh

    cd ../..

done
