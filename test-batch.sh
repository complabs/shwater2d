#!/bin/bash -l

#SBATCH -J shwater2d

#SBATCH -A edu19.summer
# BATCH --reservation=summer-2019-08-21

#SBATCH -N 1
#SBATCH -t 20:00

#SBATCH -o test-batch.out
#SBATCH -e test-batch.err

module load gcc  # tegner

#module swap PrgEnv-cray PrgEnv-gnu
#module swap PrgEnv-cray PrgEnv-intel

THREADS=64  # 48 on tegner, 64 on beskow

PROG=./shwater2d
PROG=./shwater2d_opt

echo
echo "Running $PROG"
echo

for n in `seq 0 1 $THREADS`; do 
    #$PROG $n 1000 1000 0.1; 
    #$PROG $n 1024 1024 0.1; 
    $PROG $n 2048 2048 0.1; 
done

