#!/bin/bash

#for i in {2..20}; do
#    ../bin/averagedSvp -q 1 --aqcpqc --mstart 3 --mend 3 --steps "$i" --aqcdepth 2 --inst 10 -l 0
#done

#for i in {2..50}; do
#    ../bin/averagedSvp -q 1 --aqcpqc --mstart 3 --mend 3 --steps "$i" --aqcdepth 1 --inst 10 #-l 0
#done

seq 2 50 | parallel -j4 --lb "CORE=\$((({}-2)%4)); echo Running {} on \$CORE; taskset -c \$CORE ../bin/averagedSvp -q 1 --aqcpqc --mstart 3 --mend 3 --steps {} --aqcdepth 1 --inst 10 > /dev/null 2> /dev/null"
