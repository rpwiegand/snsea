#!/bin/bash

#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rpwiegand@gmail.com
#SBATCH --job-name=rv500pop

EXPNAME=`grep "^${SLURM_ARRAY_TASK_ID}:" exp-manifest.txt | cut -d":" -f2`
STUDY=ECJ-PopPop

date
env | grep "SLURM_ARRAY"
echo "Running snsea against ${EXPNAME} ..."
module load python/python-3.8.0-gcc-9.1.0
#python3  snseaNoPopulation.py  ExperimentParams/${STUDY}/${EXPNAME}.ini > Results/${STUDY}/${EXPNAME}.raw
python3  snseaPopulation.py  ExperimentParams/${STUDY}/${EXPNAME}.ini >> Results/${STUDY}/${EXPNAME}.raw
date

