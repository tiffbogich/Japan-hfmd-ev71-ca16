#!/bin/bash

cd ../

mkdir -p results/lhs_simplex

for i in hfmd_*
do

    echo $i
    cd $i/model
    plom reduce design_lhs_simplex.json -o mle.json
    cp mle.json ../../results/lhs_simplex/${i}.json
    cd ../../

done
