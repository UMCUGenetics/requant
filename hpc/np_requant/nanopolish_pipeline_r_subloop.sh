for A in $(echo "5 10 25 50 75 90");
do
	for I in `seq 0 4`;
	do
		sbatch nanopolish_pipeline_r.sh ./vars_CpG_gluc ../refcuts_CG_test_2/p${A}_${I}.fa
	done
done
