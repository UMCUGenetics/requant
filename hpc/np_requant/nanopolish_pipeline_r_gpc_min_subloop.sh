for A in `seq 11 16`;
do
	for I in `seq 0 4`;
	do
		sbatch nanopolish_pipeline_r_gpc_min.sh ./vars_GpC_meth ../refcuts_GC_minimal/p${A}_${I}.fa
	done
done
