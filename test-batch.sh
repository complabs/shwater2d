#!/bin/bash -l

#SBATCH -J shwater2d

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-batch.out
#SBATCH -e test-batch.err

# tegner:
#module load gcc  # tegner
#THREADS=48

# beskow:
THREADS=64

PROG=./shwater2d
PROG=./shwater2d_opt

echo
echo "Running $PROG"
echo

for n in `seq 0 1 $THREADS`; do 
    #$PROG $n 1000 1000 0.1; 
    $PROG $n 2000 2000 0.1; 
done

