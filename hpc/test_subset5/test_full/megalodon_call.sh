#!/bin/bash
#SBATCH -N 1 -c 8
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=8

SUBSET=meth_CpG
path_fast5_mod=/hpc/compgen/projects/requant/requant/raw/OUD6719/fast5
path_fast5_mod_test=/hpc/compgen/projects/requant/requant/raw/OUD7493/fast5
path_fast5_can=/hpc/compgen/projects/requant/requant/raw/VER5940/fast5

file_reference_full=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/remora/asmetha/lambda_tests/lambda_phage.fa 
file_reference=${file_reference_full}
ref_base=`basename $file_reference`
file_base=${ref_base/.fa/}


# Call modified reads with trained model
mkdir -p ./${file_base}/megalodon_call/
megalodon \
        ${path_fast5_mod_test} \
        --outputs mod_mappings \
        --reference ${file_reference_full} \
        --write-mods-text \
        --processes ${num_threads} \
        --guppy-server-path ../../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-params " \
                --num_callers 4 \
                --ipc_threads 4" \
        --guppy-config dna_r9.4.1_450bps_fast.cfg \
	--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5mc CG 0 \
        --devices cuda:0 \
        --output-directory ./${file_base}/megalodon_call/mod_rmb


# Call canonical reads with trained model
megalodon \
        ${path_fast5_can} \
        --outputs mod_mappings \
        --reference ${file_reference_full} \
        --write-mods-text \
        --processes ${num_threads} \
        --guppy-server-path ../../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-params " \
                --num_callers 4 \
                --ipc_threads 4" \
        --guppy-config dna_r9.4.1_450bps_fast.cfg \
	--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5mc CG 0 \
        --devices cuda:0 \
        --output-directory ./${file_base}/megalodon_call/can_rmb

