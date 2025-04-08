#!/bin/bash
#SBATCH -N 1 -c 8
#SBATCH --mem 16G
#SBATCH --time 36:00:00

export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/

file_train_fastq_mod=./bonito_calls_chr22.fastq
file_train_bam_mod=../bonito_calls_chr22.bam
#file_reference=/hpc/compgen/projects/requant/requant/raw/cliveome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
file_reference=./GRCh38_no_alt_analysis_set_GCA_000001405.15.CpGmeth.fasta
path_np_train=./chr22_methonly/
num_threads=8
file_model_exp=./r10.4.1_400bps.cpg.9mer.model
#file_model_exp=./r10_450bps_nucleotide_9mer_expand.model

./nanopolish_precomp train \
	--train-kmers=methylated \
	--reads ${file_train_fastq_mod} \
	--bam ${file_train_bam_mod} \
	--genome ${file_reference} \
	-d ${path_np_train} \
	-t ${num_threads} \
	-i ${file_model_exp}
