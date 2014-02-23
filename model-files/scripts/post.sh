#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i/model
    pwd

    rm design.sh.* #remove cluster logs
    rm kalman kmcmc mif pmcmc simplex simul smc worker ksimplex #remove binaries

    cd C/templates
   
    make clean
    make
    make install

    cd ../../../../
   
done
