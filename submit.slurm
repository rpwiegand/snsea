#!/bin/bash

#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rpwiegand@gmail.com
#SBATCH --job-name=hammconv

EXPNAME=`grep "^${SLURM_ARRAY_TASK_ID}:" exp-manifest.txt | cut -d":" -f2`
STUDY=ECJ-NoPop

date
env | grep "SLURM_ARRAY"
echo "Running snsea against ${EXPNAME} ..."
module load python/python-3.8.0-gcc-9.1.0
python3  snseaNoPopulation.py  ExperimentParams/${ECJ-NoPop}/${EXPNAME}.ini > Results/${ECJ-NoPop}/${EXPNAME}.raw
date

