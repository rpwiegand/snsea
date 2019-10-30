#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1-00:00:00

EXPNAME=`ls -l ExperimentParams/*.ini | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -d"/" -f2 | sed s/".ini"/""/`

date
echo "Running snsea against ${EXPNAME} ..."
module load python/python-2.7.14-gcc-7.1.0
python snsea2.py  ExperimentParams/${EXPNAME}.init > Results/${EXPNAME}.out
date

