#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    plom adapt theta.json -o model/theta_lhs_simplex.json

    cp design_lhs_simplex.json model/
    cd model

    plom bootstrap theta_lhs_simplex.json design_lhs_simplex.json --engine PBS --walltime "23:30:00"
    qsub results/lhs_simplex/design.sh

    cd ../..

done
