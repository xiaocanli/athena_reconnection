#!/bin/bash
date
module load cpu cray-hdf5-parallel
module list

export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 128 -N 1 -c 2 --cpu_bind=cores ./athena -i athinput.reconnection

date
echo 'Done'
