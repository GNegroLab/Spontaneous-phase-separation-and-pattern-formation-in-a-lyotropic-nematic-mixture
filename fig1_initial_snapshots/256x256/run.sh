#!/bin/bash
#SBATCH --job-name=fig1
#SBATCH --output=fig1.out
#SBATCH --error=fig1.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=3
#SBATCH --partition=long

cc Qrho.c -o Qrho -O3 -lm
./Qrho 0.700 1.700 && ./Qrho 1.000 2.000 && ./Qrho 1.400 2.400
