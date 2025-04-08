#!/bin/bash
#SBATCH -n 8
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=8

SUBSET=meth_CpG
path_fast5_mod=../fast5_test_3/meth_CpG/
path_fast5_can=../fast5_test_3/canon/

file_reference_full=../asmetha/lambda_tests/lambda_phage.fa
# file_reference=./refcuts_CG_6_320/J02459.1_33413-35576.fa
# file_reference=./refcuts_CG_6_320/J02459.1_23100-24684.fa
file_reference=./refcuts_CG_6_128/J02459.1_23100-23669.fa
ref_base=`basename $file_reference`
file_base=${ref_base/.fa/}


# Modified reads signal mapping
mkdir -p ./${file_base}/megalodon_sigmap/
megalodon \
        ${path_fast5_mod} \
        --reference ${file_reference} \
        --output-directory ./${file_base}/megalodon_sigmap/mod \
        --outputs mappings signal_mappings \
        --num-reads 40000 \
        --guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-config dna_r9.4.1_450bps_fast.cfg \
        --devices 0 \
        --processes ${num_threads} \
        --ref-mods-all-motifs m 5mC CG 0


# Canonical reads signal mapping
megalodon \
        ${path_fast5_can} \
        --reference ${file_reference} \
        --output-directory ./${file_base}/megalodon_sigmap/can \
        --outputs mappings signal_mappings \
        --num-reads 40000 \
        --guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-config dna_r9.4.1_450bps_fast.cfg \
        --devices 0 \
        --processes ${num_threads}


# Merge all sigmapped reads
mkdir -p ./${file_base}/taiyaki_merge/
python ../taiyaki/misc/merge_mappedsignalfiles.py \
        ./${file_base}/taiyaki_merge/merged_signal_mappings.hdf5 \
        --input ./${file_base}/megalodon_sigmap/can/signal_mappings.hdf5 None \
        --input ./${file_base}/megalodon_sigmap/mod/signal_mappings.hdf5 None \
        --allow_mod_merge \
        --batch_format


# Prepare data for remora training
mkdir -p ./${file_base}/remora_prep/
remora \
        dataset prepare \
        ./${file_base}/taiyaki_merge/merged_signal_mappings.hdf5 \
        --output-remora-training-file ./${file_base}/remora_prep/remora_train_chunks.npz \
        --motif CG 0 \
        --chunk-context 50 50 \
        --kmer-context-bases 2 2 \
        --max-chunks-per-read 50 \
        --log-filename ./${file_base}/remora_prep/remora_prep_log.txt


# Actually train a model
# mkdir -p ./${file_base}/remora_train/
remora \
        model train ./${file_base}/remora_prep/remora_train_chunks.npz \
        --model ../remora/models/ConvLSTM_w_ref.py \
        --device 0 \
        --size 96 \
        --epochs 250 \
        --early-stopping 10 \
        --scheduler StepLR \
        --lr-sched-kwargs step_size 10 int \
        --lr-sched-kwargs gamma 0.5 float \
        --output-path ./${file_base}/remora_train/ \
        --device 0


# Call modified reads with trained model
mkdir -p ./${file_base}/megalodon_call/
megalodon \
        ${path_fast5_mod} \
        --outputs mod_mappings \
        --reference ${file_reference_full} \
        --remora-model ./${file_base}/remora_train/model_best.onnx \
        --write-mods-text \
        --processes ${num_threads} \
        --guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-params " \
                --num_callers 4 \
                --ipc_threads 4" \
        --guppy-config dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg \
        --devices cuda:0 \
        --output-directory ./${file_base}/megalodon_call/mod


# Call canonical reads with trained model
megalodon \
        ${path_fast5_can} \
        --outputs mod_mappings \
        --reference ${file_reference_full} \
        --remora-model ./${file_base}/remora_train/model_best.onnx \
        --write-mods-text \
        --processes ${num_threads} \
        --guppy-server-path ../ont-guppy_v5/bin/guppy_basecall_server \
        --guppy-params " \
                --num_callers 4 \
                --ipc_threads 4" \
        --guppy-config dna_r9.4.1_450bps_modbases_5mc_hac_prom.cfg \
        --devices cuda:0 \
        --output-directory ./${file_base}/megalodon_call/can