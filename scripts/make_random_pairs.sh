#!/bin/bash

cd ~/spouse_gwas/

traits="chd,height"

for i in {1..100}
do
	echo $i
	Rscript scripts/random_non_spouse_pairs.R $i $traits "chd"
	echo "Finished $i"
done 

