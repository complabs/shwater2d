#!/bin/bash -l

#SBATCH -J shwater-f90

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 30:00

#SBATCH -o test-f90-naive.out
# BATCH -e test-f90.err

module swap PrgEnv-cray PrgEnv-gnu

cd f90

THREADS=64 # on beskow

for n in `seq 1 1 $THREADS`; do 
    export OMP_NUM_THREADS=$n
    echo
    echo "Threads $n "
    srun -n 1 ./shwater2d_naive
done

