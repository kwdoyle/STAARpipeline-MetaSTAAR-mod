#!/bin/bash

saveprefix=topmed_chr_all
#saveprefix=10k_chr_all
#saveprefix=10k_10th_pcnt_telo_case_chr_all


for i in {1..22}; do
	echo topmed_chr_${i}.out >> mergeBED.list
	#echo 10k_chr_${i}.out >> mergeBED.list
	#echo 10k_10th_pcnt_telo_case_chr_${i}.out >> mergeBED.list
done

# make full bed file
echo Making full bed file
# ../plink  # originally used plink binary from FastSparseGRM, but should be able to use any plink (in theory)
plink --merge-list mergeBED.list --make-bed --out $saveprefix

# prune bed file (make lists)
echo Pruning bed file
plink  --bfile $saveprefix --indep-pairwise 50 5 0.1 --out ${saveprefix}.prunedlist

# save pruned bed
echo Saving pruned bed file
plink --bfile $saveprefix --extract ${saveprefix}.prunedlist.prune.in --make-bed --out ${saveprefix}_pruned

