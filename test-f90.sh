#!/bin/bash -l

#SBATCH -J shwater-f90

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-f90.out
# BATCH -e test-f90.err

cd f90

echo
echo "Running f90"
echo

OMP_NUM_THREADS=32 ./shwater2d

