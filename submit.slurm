#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1-00:00:00

N=64
RHOMIN=6
MAXGEN=100000

python snsea.py  --rhoMin=$RHOMIN  --n=$N  --maxGen=$MAXGEN  --archiveDistancesFile=$RHOMIN-$N.out
