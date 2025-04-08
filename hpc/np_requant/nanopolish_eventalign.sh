#!/bin/bash
#SBATCH -n 4
#SBATCH --mem 16G
#SBATCH --time 24:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=4

file_reference_full=../../asmetha/lambda_tests/lambda_phage.fa
file_fastq_mod=./guppy_basecall_test_4/meth_CpG/fastq_runid_459de6c83c60c638bb74f3e431d101f3942d1fac.fastq
file_bam_mod=./minimap2/lambda_phage_meth_CpG_sort.bam
file_fastq_canon=./guppy_basecall_test_4/canon/fastq_runid_248f8352d7c1a9aaaeb9f5207ceb43108a8cd9cb.fastq
file_bam_canon=./minimap2/lambda_phage_canon_sort.bam


file_eventalign=./test_eventalign/${file_base} # Base file for output data

mkdir -p ./test_eventalign/

#./nanopolish_og eventalign \
#	--reads   ${file_fastq_mod} \
#	--bam     ${file_bam_mod} \
#	--genome  ${file_reference_full} \
#	--scale-events > ${file_eventalign}/mod.eventalign.txt

./nanopolish_og eventalign \
	--reads   ${file_fastq_canon} \
	--bam     ${file_bam_canon} \
	--genome  ${file_reference_full} \
	--scale-events > ${file_eventalign}/canon.eventalign.txt
