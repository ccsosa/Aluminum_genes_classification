#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=30
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load lang/python/3.9.7
python /users/ccsosaa/pecanpy/PULEARN_ALUMINUM.py