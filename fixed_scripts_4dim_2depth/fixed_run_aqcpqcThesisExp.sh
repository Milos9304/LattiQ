#!/bin/bash

#for i in {2..20}; do
#    ../bin/averagedSvp -q 1 --aqcpqc --mstart 3 --mend 3 --steps "$i" --aqcdepth 2 --inst 10 -l 0
#done

#for i in {2..50}; do
#    ../bin/averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps "$i" --aqcdepth 2 --inst 10 -l 0
#done

#seq 2 50 | parallel -j4 --ungroup --lb taskset -c {#}%4 ../bin/averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps {} --aqcdepth 2 --inst 10 -l 0
seq 2 50 | parallel -j4 --lb "CORE=\$((({}-2)%4)); echo Running {} on \$CORE; taskset -c \$CORE ../bin/averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps {} --aqcdepth 2 --inst 10 -l 0 > /dev/null 2> /dev/null"

