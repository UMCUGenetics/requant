for A in `seq 1 16`;
do
	for I in `seq 0 4`;
	do
		sbatch nanopolish_pipeline_r_min.sh ./vars_CpG_meth ../refcuts_CG_minimal_2/p${A}_${I}.fa
#		sbatch nanopolish_pipeline_r_min.sh ./vars_CpG_gluc ../refcuts_CG_minimal/p${A}_${I}.fa
	done
done
