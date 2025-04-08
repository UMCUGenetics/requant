for A in $(echo "5 10 25 50 75 90");
do
	for I in `seq 0 4`;
	do
		sbatch nanopolish_pipeline_r_gluc_on_meth.sh ./vars_CpG_gluc_on_meth ../refcuts_CG_test_2/p${A}_${I}.fa
	done
done
