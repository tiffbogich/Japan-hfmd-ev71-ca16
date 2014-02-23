#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i
    plom build
    cp ../results/lhs_simplex/$i.json model/theta.json

    cd model
    fit replicate theta.json "smc ode" 1 --run

    cd ../
    plom commit model/design.json

    cd ../

done
