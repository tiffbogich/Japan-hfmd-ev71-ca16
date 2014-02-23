#!/bin/bash

cd ../

for i in hfmd_*
do
    echo $i model:

    cd $i
    plom build --local
    cd ../

done
