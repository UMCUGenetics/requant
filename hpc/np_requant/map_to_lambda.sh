conda activate poretools

FASTQ_IN=$1
PATH_OUT=$2

minimap2 -a -x map-ont -t 2 \
    -o ${PATH_OUT}.sam \
    /hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/remora/asmetha/lambda_tests/lambda_phage.fa \
    ${FASTQ_IN}

samtools_1.17 view -b ${PATH_OUT}.sam > ${PATH_OUT}.bam
samtools_1.17 sort ${PATH_OUT}.bam > ${PATH_OUT}_sort.bam
samtools_1.17 index ${PATH_OUT}_sort.bam

rm ${PATH_OUT}.sam ${PATH_OUT}.bam