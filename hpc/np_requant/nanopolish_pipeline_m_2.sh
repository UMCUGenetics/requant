for A in $(echo "5 10 25 50 75 90");
do
	for I in `seq 1 4`;
	do
		sbatch nanopolish_pipeline_m.sh ../refcuts_CG_test_2/p${A}_${I}.fa
	done
done
