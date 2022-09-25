cmake .
make
for SIMDLEN in 1 4 8 16
do
	for data in 'afiro' '25fv47' 'degen2' 'Recon1'
	do
		./crhmc_sample_sparse $data $SIMDLEN 80000 20000
	done
done
