#!/bin/bash
#SBATCH -n 4
#SBATCH --mem 16G
#SBATCH --time 24:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=4

file_reference_full=../../asmetha/lambda_tests/lambda_phage.fa
file_fastq_mod=./guppy_basecall/fastq_runid_459de6c83c60c638bb74f3e431d101f3942d1fac.fastq
file_bam_mod=./minimap2/lambda_phage_pass_sort.bam
file_fastq_canon=./guppy_basecall_canon/fastq_runid_248f8352d7c1a9aaaeb9f5207ceb43108a8cd9cb.fastq
file_bam_canon=./minimap2_canon/lambda_phage_canon_pass_sort.bam


# file_reference=../refcuts_CG_6_320/J02459.1_33413-35576.fa
# file_reference=../refcuts_CG_6_320/J02459.1_23100-24684.fa
file_reference=../refcuts_CG_6_128/J02459.1_23100-23669.fa
ref_base=`basename $file_reference`
file_base=${ref_base/.fa/}

file_bam=./${file_base}/minimap2/${file_base}
path_np_train=./${file_base}/nanopolish_train/
file_rq=./${file_base}/requant/${file_base}

# Setup model file names and variables
file_model_exp=./r9.4_450bps.nucleotide.6mer.expand_cpg.model # Expanded nanopolish model
file_model_dir=${path_np_train}/r9.4_450bps.cpg.6mer.template.round4.model # Model as trained (directly)
file_model_rep=${file_rq}.replaced.model # Model after replacing all values with imputed
file_model_add=${file_rq}.added.model # Model after filling missing values with imputed

file_methcalls=./${file_base}/nanopolish_methcall/${file_base} # Base file for output data


mkdir -p ./${file_base}/minimap2/ ./${file_base}/nanopolish_train/ ./${file_base}/requant ./${file_base}/nanopolish_methcall/

### Region specific steps ###

## Map modified reads to reference genome section ##
minimap2 -a \
    -x map-ont \
    -o ${file_bam}.sam \
    ${file_reference} \
    ${file_fastq_mod}

# Sam to bam, sort and index
samtools view -b ${file_bam}.sam > ${file_bam}.bam
samtools sort ${file_bam}.bam > ${file_bam}.sort.bam
samtools index ${file_bam}.sort.bam
rm ${file_bam}.sam ${file_bam}.bam
file_bam=${file_bam}.sort.bam


## Train nanopolish model ##
# Methylate all motif matches in reference sequence
python methylate_reference.py \
    --recognition cpg \
    ${file_reference} \
    > ${file_reference/.fa/_meth_cpg.fa}

# Train nanopolish model
./nanopolish_og train \
    --reads ${file_fastq_mod} \
    --bam ${file_bam} \
    --genome ${file_reference/.fa/_meth_cpg.fa} \
    -d ${path_np_train} \
    -t ${num_threads} \
    -i ${file_model_exp}


## Impute kmers ##
python requant_impute.py \
    ${path_np_train}/r9.4_450bps.cpg.6mer.template.round4.model \
    ${file_reference/.fa/.log} \
    ${file_rq}


### Call methylation ###

## Replaced model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_rep} \
	--reads   ${file_fastq_mod} \
	--bam     ${file_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.replaced.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_rep} \
	--reads   ${file_fastq_canon} \
	--bam     ${file_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.replaced.canon.methcalls


## Replaced model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_add} \
	--reads   ${file_fastq_mod} \
	--bam     ${file_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.added.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_add} \
	--reads   ${file_fastq_canon} \
	--bam     ${file_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.added.canon.methcalls
    

## Not imputed model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_dir} \
	--reads   ${file_fastq_mod} \
	--bam     ${file_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.direct.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_dir} \
	--reads   ${file_fastq_canon} \
	--bam     ${file_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.direct.canon.methcalls




# # Reference function for basecalling data, only needed once per dataset 
# # (i.e. one run for modified data and one for canonical data).
# guppy_basecall(){
#     local path_fast5=$1
#     local path_fastq=$2

#     # Basecalling modified data
#     ../../ont-guppy_v5/bin/guppy_basecaller \
#         -i ${path_fast5} \
#         -s ${path_fastq} \
#         -c dna_r9.4.1_450bps_modbases_5mc_hac.cfg \
#         -x cuda:0

#     # Merge reads to single fastq
#     cat ${path_fastq}/*/*.fastq > ${path_fastq}/merged.fastq

#     # Index the reads for nanopolish
#     ./nanopolish_og index \
#         -d ${path_fast5} \
#         ${path_fastq}/merged.fastq
# }


### Not specific to region ###

# Expand model
#python expand_model_alphabet.py --alphabet cpg /hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model > r9.4_450bps.nucleotide.6mer.expand_cpg.model

# Basecalling all data for training and testing
#guppy_basecall ../../fast5_test_3/meth_CpG/ ./guppy_basecall_meth_CpG/ # Canonical
#guppy_basecall ../../fast5_test_3/canon/ ./guppy_basecall_canon/ # Modified

# Map reads to full ref genome for methylation calling later
#map_to_ref # Canonical
#map_to_ref # Modified