#! /bin/bash

fit slice theta.json "./smc deter -J 1 -h -r" 100 --run --par r0_1:all,r0_2:all,v:all,q:all,e:all,d:all,z:all,sigma:all,iota_1:all,iota_2:all,SS:all,IS:all,SI:all,SR:all,RS:all,rep:EV71_JAP__IASR__inc,rep:CA16_JAP__IASR__inc
