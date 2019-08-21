#!/bin/bash -l

#SBATCH -J shwater-ref

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-ref.out
#SBATCH -e test-ref.err

# The variable PROBLEMSIZE is globally defined
. problem-size

cd orig

THREADS=64 # on beskow

for n in `seq 1 1 $THREADS`; do 
    export OMP_NUM_THREADS=$n
    echo
    echo "Threads $n "
    srun -n 1 ./shwater2d_naive
done

