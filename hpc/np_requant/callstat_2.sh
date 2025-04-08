#!/bin/bash
#SBATCH -N 1 -c 1
#SBATCH --mem 8G
#SBATCH --time 4:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

python callstat_2.py $1 $2 > $1.stat_2.txt
