#!/bin/bash
#SBATCH -n 4
#SBATCH --mem 64G
#SBATCH --time 24:00:00
##SBATCH --partition=gpu
##SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

SUBSET=can_CpG
mkdir -p ./remora_prep/${SUBSET}/

remora \
	dataset prepare \
	./taiyaki_merge/${SUBSET}/merged_signal_mappings.hdf5 \
	--output-remora-training-file ./remora_prep/${SUBSET}/remora_train_chunks.npz \
	--motif CG 0 \
	--chunk-context 50 50 \
	--kmer-context-bases 2 2 \
	--max-chunks-per-read 50 \
	--log-filename remora_prep/${SUBSET}/remora_prep_log.txt

