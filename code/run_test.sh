#!/bin/bash
date
module swap craype-haswell craype-mic-knl
module load cray-hdf5-parallel
module list

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export MPICH_MAX_THREAD_SAFETY=multiple

# export PMI_MMAP_SYNC_WAIT_TIME=400
time srun -n 256 -N 4 -c 4 --cpu_bind=cores ./athena -i athinput.reconnection

date
echo 'Done'
