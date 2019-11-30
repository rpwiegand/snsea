#!/bin/bash

#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=1-00:00:00

EXPNAME=hamming-r3-n16

date
echo "Running snsea against ${EXPNAME} starting on trial number ${SLURM_ARRAY_TASK_ID} ..."
module load python/python-2.7.14-ic-2017.1.043
python snseaPopulation.py  ExperimentParams/${EXPNAME}.ini --startTrialNum=${SLURM_ARRAY_TASK_ID} >> Results/${EXPNAME}.raw
date

