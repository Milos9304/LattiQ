# LattiQ
Variational quantum algorithms experimental framework for lattice problems based on [FastVQA](https://github.com/Milos9304/FastVQA). The experiments are reported in [Variational quantum solutions to the Shortest Vector Problem](https://arxiv.org/abs/2202.06757?fbclid=IwAR1BMNjJZ2KCKjE3vckYgiRg4V5hE-aDMDIpO9CqDwwM8tAN7tD9PW1QojU).

## TODO
- Allow for mappings that result in multiple zero reference states.

## Features
- loads matrices in experiments folder
- encodes the Shortest Vector Problem (SVP) to an Ising spin Hamiltonian
- runs VQE/QAOA algorithm to find the shortest lattice vector and reports the results

## Installation
```
mkdir build && cd build
cmake .. && make
make
```

## Usage
```
./bin/LattiQ -e --vqe -r 25 -x 0.175
```
The **-e** option says that the qary_25_50 experiment should be run, **vqe**/**qaoa** specifies an algorithm to run, **-r** specifies the rank for the experiment as described in the paper and **-x** is CVaR's alpha value.
