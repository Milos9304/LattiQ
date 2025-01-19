#!/bin/bash

#for i in {2..20}; do
#    ./averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps "$i" --aqcpqcdepth 1 --inst 10 
#done

for i in {2..20}; do
    ./averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps "$i" --aqcdepth 2 --inst 10
done
