#!/bin/bash
date
module load cpu cray-hdf5-parallel
module list

export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

time srun -n 256 -N 4 -c 4 --cpu_bind=cores ./athena -i $1

date
echo 'Done'
