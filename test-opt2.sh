#!/bin/bash -l

#SBATCH -J shwater-opt

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-opt2.out
# BATCH -e test-opt2.err

# tegner:
#module load gcc  # tegner
#THREADS=48

PROG=./shwater2d_opt

echo "C++, optimized, 8 threads"

for n in `seq 500 100 2000`; do 
    $PROG 8 $n $n 0.1
done

