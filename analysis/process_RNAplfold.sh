#!/bin/bash

mkdir RNAplfold_results
cat target_sequence.75.csv | cut -d ',' -f 1,2 | { read header; while IFS="," read p seq
do
	echo $p >> RNAplfold_results.txt

	echo \>${p} > ./$p.fa
	echo $seq >> ./$p.fa
	cat ./$p.fa | RNAplfold -L 40 -W 80 -u 50
	cp ./${p}_lunp ./RNAplfold_results/
	cp ./${p}_dp.ps ./RNAplfold_results/
	
	rm $p.fa
	rm ./${p}_lunp
	rm ./${p}_dp.ps
	echo $p Done!
done; }