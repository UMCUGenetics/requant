#!/bin/bash
#SBATCH -N 1 -c 4
#SBATCH --mem 4G
#SBATCH --time 1:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/

num_threads=4

A=50
I=0
file_reference=../refcuts_GC_test_2/p${A}_${I}.fa
sample_mod_train=OUD6724
file_train_fastq_mod=../../../raw/${sample_mod_train}/${sample_mod_train}_all.fastq
file_train_bam_mod=../../../raw/${sample_mod_train}/${sample_mod_train}_all_sort.bam
file_model_exp=./r9.4_450bps.nucleotide.6mer.expand_gpc.model

path_np_train=./nanopolish_train_debug


./nanopolish_og train \
    --train-kmers=methylated \
    --reads ${file_train_fastq_mod} \
    --bam ${file_train_bam_mod} \
    --genome ${file_reference} \
    -d ${path_np_train} \
    -t ${num_threads} \
    -i ${file_model_exp}
