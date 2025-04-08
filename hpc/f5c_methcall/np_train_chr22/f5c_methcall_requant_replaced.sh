#!/bin/bash
#SBATCH -n 8
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --partition=gpu
#SBATCH --gpus-per-node=1

source /hpc/compgen/users/rstraver/miniconda3/bin/activate poretools

# Fix missing VBZ plugin
export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/

PATH_CLIVE=/hpc/compgen/projects/requant/requant/raw/cliveome/

../../f5c-v1.3/f5c_x86_64_linux_cuda call-methylation \
    --iop 4 \
    -b ../bonito_calls_chr21.bam \
    -g ${PATH_CLIVE}/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -r ./bonito_calls_chr21.fastq \
    --kmer-model ../tnmm_r10_450bps.nucleotide.9mer.template.round9.model \
    --meth-model ./chr22_methonly/rq_test.replaced.model \
        > ./chr22_methonly/PAM63167_result_requant_replaced.tsv
