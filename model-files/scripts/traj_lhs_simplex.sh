#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/lhs_simplex/$i.json model/mle.json

    cd model
    plom pipe mle.json | ./smc ode -J 100 --traj

    cd ../..

done
