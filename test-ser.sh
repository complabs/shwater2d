#!/bin/bash -l

#SBATCH -J shwater-ser

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-ser.out
# BATCH -e test-ser.err

# tegner:
#module load gcc  # tegner
#THREADS=48

# beskow:
THREADS=64

# The variable PROBLEMSIZE is globally defined
. problem-size

PROG=./shwater2d_ser

echo
echo "Running $PROG, size: $PROBLEMSIZE"
echo

$PROG 0 $PROBLEMSIZE

