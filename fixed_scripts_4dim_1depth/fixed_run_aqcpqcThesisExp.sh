#!/bin/bash

#for i in {2..20}; do
#    ../bin/averagedSvp -q 1 --aqcpqc --mstart 3 --mend 3 --steps "$i" --aqcdepth 2 --inst 10 -l 0
#done

for i in {15..50}; do
    ../bin/averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps "$i" --aqcdepth 1 --inst 10 -l 0
done
