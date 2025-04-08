#!/bin/bash
#SBATCH -n 10
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

SUBSET=meth_CpG
mkdir -p ./megalodon_sigmap/

megalodon \
	../fast5_test_3/${SUBSET}/ \
	--reference ./refcuts_CG_6_320/J02459.1_33413-35576.fa \
	--output-directory ./megalodon_sigmap/${SUBSET}/ \
	--outputs mappings signal_mappings \
	--num-reads 40000 \
	--guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
	--guppy-config dna_r9.4.1_450bps_fast.cfg \
	--devices 0 \
	--processes 10 \
	--ref-mods-all-motifs m 5mC CG 0
