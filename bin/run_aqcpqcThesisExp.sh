#!/bin/bash

for i in {5..20}; do
    ./averagedSvp -q 1 --aqcpqc --mstart 4 --mend 4 --steps "$i" --aqcpqcdepth 1 --inst 2
done
