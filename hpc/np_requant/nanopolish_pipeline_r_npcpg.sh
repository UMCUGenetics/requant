#!/bin/bash
#SBATCH -N 1 -n 4
#SBATCH --mem 4G
#SBATCH --time 36:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=4

# Load sample specific variables
source ./vars_CpG_meth

file_reference=./event_align/lambda_phage_methylated_cpg.fa
ref_base=`basename $file_reference`
file_base=${ref_base/.fa/}
path_base="np_cpg"
path_base=`echo "${path_output}/${path_base}/${file_base}"`

path_np_train=./${path_base}/nanopolish_train/
file_rq=./${path_base}/requant/${file_base}

# Setup model file names and variables
file_model_exp=./r9.4_450bps.nucleotide.6mer.expand_cpg.model # Expanded nanopolish model
file_model_dir=./r9.4_450bps.cpg.6mer.template.model # Model as trained (directly)
file_model_rep=${file_rq}.replaced.model # Model after replacing all values with imputed
file_model_add=${file_rq}.added.model # Model after filling missing values with imputed

file_methcalls=${path_base}/nanopolish_methcall/${file_base} # Base file for output data


mkdir -p ${path_base}/nanopolish_train/ ${path_base}/requant ${path_base}/nanopolish_methcall/

### Region specific steps ###

# # Train nanopolish model 
# # Use already mapped full data and already M added reference file
# ./nanopolish_og train \
#    --train-kmers=methylated \
#    --reads ${file_train_fastq_mod} \
#    --bam ${file_train_bam_mod} \
#    --genome ${file_reference} \
#    -d ${path_np_train} \
#    -t ${num_threads} \
#    -i ${file_model_exp}


## Impute kmers ##
python requant_impute.py \
    ${file_model_dir} \
    - \
    ${file_rq} \
		>	${file_rq}.rules


# ### Call methylation ###

## Replaced model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_rep} \
	--reads   ${file_test_fastq_mod} \
	--bam     ${file_test_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.replaced.mod.methcalls
bzip2 ${file_methcalls}.replaced.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_rep} \
	--reads   ${file_test_fastq_canon} \
	--bam     ${file_test_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.replaced.canon.methcalls
bzip2 ${file_methcalls}.replaced.canon.methcalls


## Added model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_add} \
	--reads   ${file_test_fastq_mod} \
	--bam     ${file_test_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.added.mod.methcalls
bzip2 ${file_methcalls}.added.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_add} \
	--reads   ${file_test_fastq_canon} \
	--bam     ${file_test_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.added.canon.methcalls
bzip2 ${file_methcalls}.added.canon.methcalls
    

## Not imputed model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_dir} \
	--reads   ${file_test_fastq_mod} \
	--bam     ${file_test_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.direct.mod.methcalls
bzip2 ${file_methcalls}.direct.mod.methcalls
# Canonical
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_dir} \
	--reads   ${file_test_fastq_canon} \
	--bam     ${file_test_bam_canon} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.direct.canon.methcalls
bzip2 ${file_methcalls}.direct.canon.methcalls
