#!/bin/bash -l

#SBATCH -J shwater-opt

#SBATCH -A edu19.summer
#SBATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-opt.out
# BATCH -e test-opt.err

# tegner:
#module load gcc  # tegner
#THREADS=48

# beskow:
THREADS=64

# The variable PROBLEMSIZE is globally defined
. problem-size

PROG=./shwater2d_opt

echo
echo "Running $PROG, size: $PROBLEMSIZE"
echo

for n in `seq 0 1 $THREADS`; do 
    $PROG $n $PROBLEMSIZE
done

