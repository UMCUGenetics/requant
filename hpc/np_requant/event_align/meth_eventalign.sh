set -e
#conda activate poretools


# FOLD_NAME=OUD7495
# FOLD_PASS=pass
# FILE_NAME=AID494_${FOLD_PASS}_6249157e_34
# TYPE_NAME=gpc_${FILE_NAME}

#../../../../raw/OUD7495/fastq_fail/AJG810_fail_c4241646_24.fastq.gz
FOLD_NAME=OUD7495
FOLD_PASS=fail
FILE_NAME=AJG810_${FOLD_PASS}_c4241646_$1
TYPE_NAME=glu_${FILE_NAME}

export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/
# FILE_REFERENCE=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/remora/asmetha/lambda_tests/lambda_phage.fa
FILE_REFERENCE=./lambda_phage_methylated_cpg.fa

# echo make symlink
# ln -s ../../../../raw/${FOLD_NAME}/fastq_${FOLD_PASS}/${FILE_NAME}.fastq.gz ${TYPE_NAME}.fastq.gz
## mv ${FILE_NAME}.fastq.gz ${TYPE_NAME}.fastq.gz

#echo index reads
#../nanopolish_og index -d ../../../../raw/${FOLD_NAME}/fast5/ ${TYPE_NAME}.fastq.gz

#echo map to lamba
#./map_to_lambda.sh ./${TYPE_NAME}.fastq.gz ./${TYPE_NAME}

echo event align
# ./eventalign.sh ./${TYPE_NAME}.fastq.gz ./${TYPE_NAME}.bam > ./${TYPE_NAME}_eventalign.tsv
../nanopolish_og \
        eventalign \
        -r ./${TYPE_NAME}.fastq.gz \
        -b ./${TYPE_NAME}_sort.bam \
        -g ${FILE_REFERENCE} \
        -t 4 \
        --scale-events \
        --samples \
        --signal-index \
        > ./${TYPE_NAME}_eventalign_mcpg.tsv

echo bzip2 compress
bzip2 ./${TYPE_NAME}_eventalign_mcpg.tsv




