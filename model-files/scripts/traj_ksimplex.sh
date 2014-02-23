#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/ksimplex/$i.json model/mle.json

    cd model
    plom pipe mle.json | ./kalman --no_dem_sto --traj

    cd ../..

done
