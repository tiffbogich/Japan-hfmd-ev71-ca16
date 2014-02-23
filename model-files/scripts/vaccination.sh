#!/bin/bash


for i in mle/*.json
do

    myfile=$(basename "$i")
    extension="${myfile##*.}"
    filename="${myfile%.*}"

    cd ${filename}_v/
    pwd

    pmbuilder context.json process.json link.json -o model && fit adapt theta.json context.json process.json link.json > model/theta.json
    cp design.json model/

    cd model
    fit merge theta.json ../../mle/${myfile} -o theta.json

    fit bootstrap theta.json design.json --run


    cd ../../

done
