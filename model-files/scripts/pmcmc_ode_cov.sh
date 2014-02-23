#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/lhs_simplex/$i.json model/theta_pmcmc_ode_cov.json
    cp design_pmcmc_ode_cov.json model/

    cd model
    plom bootstrap theta_pmcmc_ode_cov.json design_pmcmc_ode_cov.json --engine PBS --walltime 72:00:00
    qsub results/pmcmc_ode_cov/design.sh

    cd ../..

done
