#!/bin/bash
#
#SBATCH -q regular
#SBATCH -N 4
#SBATCH -t 24:00:00
#SBATCH -C knl,quad,cache
#SBATCH -o athena%j.out
#SBATCH -e athena%j.err
#SBATCH -J reconnection
#SBATCH -A m4054
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=phyxiaolee@gmail.com
#SBATCH -L SCRATCH,project

##### These are shell commands
date
module swap craype-haswell craype-mic-knl
# module unload craype-hugepages2M
module load cray-hdf5-parallel
# module load lustre-default
module list

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export MPICH_MAX_THREAD_SAFETY=multiple

#Quad Cache:
# export PMI_MMAP_SYNC_WAIT_TIME=400
time srun -n 256 -N 4 -c 4 --cpu_bind=cores ./athena -i athinput.reconnection

date
echo 'Done'
