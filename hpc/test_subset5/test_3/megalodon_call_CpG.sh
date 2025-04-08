#!/bin/bash
#SBATCH -n 10
#SBATCH --mem 128G
#SBATCH --time 144:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

SUBSET=can_CpG
mkdir -p ./megalodon_call/${SUBSET}/

megalodon \
	/hpc/compgen/projects/asmethylation/raw/OUD7898/fast5_pass/ \
	--outputs mod_mappings \
	--reference /hpc/compgen/users/cvermeulen/refgenomes/nohap_varaan/hg19.fa \
	--remora-model ./remora_train/${SUBSET}/model_best.onnx \
	--write-mods-text \
	--processes 10 \
	--guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
	--guppy-params " \
		--num_callers 5 \
		--ipc_threads 6" \
	--guppy-config dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg \
	--devices cuda:0 \
	--output-directory ./megalodon_call/${SUBSET}/prom_CpG_pass

