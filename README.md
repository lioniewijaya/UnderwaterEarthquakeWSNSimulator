# WSNSimulatorForUnderwaterEarthquakes

This project is developed for FIT3143 Parallel Computing - S2 2022 Monash University. It is simulation of a distributed wireless sensor network (WSN) of a set of interconnected seismic sensors to detect underwater earthquakes. The simulation model is designed and implemented with parallel computing architecture using Message Passing Interface (MPI) and POSIX in C.

## Collaborators

-   Editha Karina Hermawan
-   Lionie Annabella Wijaya

## How to Run

_To compile:_<br>
`      make all
     `
<br>or<br>
`      mpicc wsn.c utils.c alert.c balloon.c node.c -lm -lpthread -o wsn
     `
<br>
<br>_To run via command line argument in VM:_<br>
version 1 - fixed iterations
<br>
`      mpirun --oversubscribe -np [number of process] wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold] [number of iterations]
     `
<br>
version 2 - programs stops during runtime after user set a sentinel value of 2 in sentinel.txt
<br>
`      mpirun --oversubscribe -np [number of process] wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold]
     `
<br>
<br>_To run via CAAS:_<br>
`      srun wsn [grid row] [grid column] [difference in magnitude threshold] [difference in distance threshold] [number of iterations]
     `

## Documents

Refer to [report](report.pdf).
