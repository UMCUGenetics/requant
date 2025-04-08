#!/bin/bash
#SBATCH -n 2
#SBATCH --mem 16G
#SBATCH --time 24:00:00
##SBATCH --partition=gpu
##SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

mkdir -p ./taiyaki_merge/can_CpG/

python ../taiyaki/misc/merge_mappedsignalfiles.py \
	./taiyaki_merge/can_CpG/merged_signal_mappings.hdf5 \
	--input ./megalodon_sigmap/canon/signal_mappings.hdf5 None \
	--input ./megalodon_sigmap/meth_CpG/signal_mappings.hdf5 None \
	--allow_mod_merge \
	--batch_format
