#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=18
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load lang/r/4.2.1  
Rscript /users/ccsosaa/pecanpy/AdaSampling_UNWEIGHTED_JJG.R