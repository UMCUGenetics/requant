export HDF5_PLUGIN_PATH=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/nanopolish/hdf5_plugin/

file_reference_full=/hpc/compgen/projects/asmethylation/asmethylation/analysis/rstraver/remora/asmetha/lambda_tests/lambda_phage.fa

sample_mod_train=OUD6719
sample_mod_test=OUD7493
sample_can_test=VER5940


export file_train_fastq_mod=../../../raw/${sample_mod_train}/${sample_mod_train}_all.fastq
export file_test_fastq_mod=../../../raw/${sample_mod_test}/${sample_mod_test}_all.fastq
export file_test_fastq_canon=../../../raw/${sample_can_test}/${sample_can_test}_all.fastq

export file_train_bam_mod=../../../raw/${sample_mod_train}/${sample_mod_train}_all_sort.bam
export file_test_bam_mod=../../../raw/${sample_mod_test}/${sample_mod_test}_all_sort.bam
export file_test_bam_canon=../../../raw/${sample_can_test}/${sample_can_test}_all_sort.bam


../nanopolish_og \
	eventalign \
	-r $1 \
	-b $2 \
	-g ${file_reference_full} \
	-t 4 \
	--scale-events \
	--samples \
	--signal-index
