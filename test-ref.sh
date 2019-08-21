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

echo
echo "Running shwater2d $PROBLEMSIZE"
echo

./shwater2d $PROBLEMSIZE

echo
echo "Running 1 thread shwater2d_naive $PROBLEMSIZE"
echo

OMP_NUM_THREADS=1 ./shwater2d_naive $PROBLEMSIZE

echo
echo "Running 16 thread shwater2d_naive $PROBLEMSIZE"
echo

OMP_NUM_THREADS=16 ./shwater2d_naive $PROBLEMSIZE

echo
echo "Running 32 shwater2d_naive $PROBLEMSIZE"
echo

OMP_NUM_THREADS=32 ./shwater2d_naive $PROBLEMSIZE

echo
echo "Running def #threads shwater2d_naive $PROBLEMSIZE"
echo

./shwater2d_naive $PROBLEMSIZE
