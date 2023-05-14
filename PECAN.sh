#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o %j.out
##SBATCH -e %j.err

#SBATCH --cpus-per-task=1
#SBATCH --partition=FULL
##SBATCH --nodelist=node15

module load  lang/python/3.9.7
pecanpy --input /users/ccsosaa/pecanpy/BIG_COMP_W.edge --output /users/ccsosaa/pecanpy/BIG_COMP_W.emb --mode DenseOTF --extend