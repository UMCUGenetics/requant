#!/bin/bash
#SBATCH -n 2
#SBATCH --mem 16G
#SBATCH --time 24:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

SUBSET=can_CpG
mkdir -p ./remora_train/

remora \
	model train ./remora_prep/${SUBSET}/remora_train_chunks.npz \
	--model ../remora/models/ConvLSTM_w_ref.py \
	--device 0 \
	--size 96 \
	--epochs 250 \
	--early-stopping 10 \
	--scheduler StepLR \
	--lr-sched-kwargs step_size 10 int \
	--lr-sched-kwargs gamma 0.5 float \
	--output-path ./remora_train/${SUBSET}
