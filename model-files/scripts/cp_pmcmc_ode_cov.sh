#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/lhs_simplex/$i.json model/theta_pmcmc_ode.json
    cp design_pmcmc_ode.json model/

    cd model

    fit bootstrap theta_pmcmc_ode.json design_pmcmc_ode.json --queue PBS,71:00:00
    mv run_designs.sh run_pmcmc_ode.sh

    cp results/pmcmc_ode_cov/*.csv results/pmcmc_ode/

    cd ../..

done
