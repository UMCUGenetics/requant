#!/bin/bash
#SBATCH -N 1 -c 4
#SBATCH --mem 4G
#SBATCH --time 36:00:00

source /hpc/compgen/users/rstraver/miniconda3/bin/activate remora

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
num_threads=4

# Load sample specific variables
source $1

file_reference=$2
ref_base=`basename $file_reference`
file_base=${ref_base/.fa/}
path_base_tmp=`echo ${file_reference} | rev | cut -d/ -f2 | rev`
path_base=`echo "${path_output}/${path_base_tmp}/${file_base}"`
path_base_gluc=`echo "./output_CpG_gluc/${path_base_tmp}/${file_base}"`
path_np_train=./${path_base_gluc}/nanopolish_train/
file_rq=./${path_base_gluc}/requant/${file_base}

# Setup model file names and variables
# file_model_exp=./r9.4_450bps.nucleotide.6mer.expand_cpg.model # Expanded nanopolish model
file_model_dir=${path_np_train}/r9.4_450bps.cpg.6mer.template.round4.model # Model as trained (directly)
file_model_rep=${file_rq}.replaced.model # Model after replacing all values with imputed
file_model_add=${file_rq}.added.model # Model after filling missing values with imputed

file_methcalls=${path_base}/nanopolish_methcall/${file_base} # Base file for output data


mkdir -p ${path_base}/nanopolish_train/ ${path_base}/requant ${path_base}/nanopolish_methcall/

### Call methylation ###

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


## Replaced model
# Modified
./nanopolish_og call-methylation \
	--threads ${num_threads} \
	--model   ${file_model_add} \
	--reads   ${file_test_fastq_mod} \
	--bam     ${file_test_bam_mod} \
	--genome  ${file_reference_full} \
	    >     ${file_methcalls}.added.mod.methcalls
bzip2 ${file_methcalls}.added.mod.methcalls


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

