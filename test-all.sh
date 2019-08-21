#!/bin/bash

echo
echo "make clean ---------------------"
echo

make clean

echo
echo "make ---------------------"
echo

make

echo

sbatch test-ref.sh
sbatch test-ser.sh
sbatch test-par.sh
sbatch test-opt.sh

