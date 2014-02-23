#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i

    cp ../results/lhs_simplex/$i.json model/theta_ksimplex.json
    cp design_ksimplex.json model/

    cd model
    plom bootstrap theta_ksimplex.json design_ksimplex.json --engine PBS --walltime "72:00:00"
    qsub results/ksimplex/design.sh

    cd ../..

done
