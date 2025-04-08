#!/bin/bash
#SBATCH -N 1 -c 1
#SBATCH --mem 4G
#SBATCH --time 2:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

python callstat.py $1 $2 > $1.stat.txt
