for A in $(echo "5 10 25 50 75 90");
do
	for I in `seq 0 4`;
	do
		sbatch nanopolish_pipeline_r_gpc.sh ./vars_GpC_meth ../refcuts_GC_test_2/p${A}_${I}.fa
	done
done
