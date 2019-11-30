#!/bin/bash

#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rpwiegand@gmail.com

#EXPNAME=`ls -l ExperimentParams/*.ini | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -d"/" -f2 | sed s/".ini"/""/`
EXPNAME=`grep "^${SLURM_ARRAY_TASK_ID}:" exp-manifest.txt | cut -d":" -f2`

date
echo "Running snsea against ${EXPNAME} ..."
module load python/python-2.7.14-ic-2017.1.043
python snseaPopuluation.py  ExperimentParams/${EXPNAME}.ini > Results/${EXPNAME}.raw
date

