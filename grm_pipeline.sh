#flprefix=${flprefix}
#flprefix=10k_chr_all
#flprefix=10k_10th_pcnt_telo_case_chr_all
flprefix=topmed_chr_all
#basedir=~/noncoding_telo/STAAR/TOPMed_Peddy_Euro/GRM/
#basedir=~/noncoding_telo/STAAR/10k_Cohort/GRM/
#basedir=~/noncoding_telo/STAAR/10k_Cohort/Telomere_Analysis_new/GRM/
basedir=~/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/GRM/
savedir=${basedir}/output/

mkdir $savedir

if [[ $1 == "king" ]]; then
	echo running king
	${basedir}/king -b ${basedir}/beds/${flprefix}.bed --ibdseg --degree 4 --cpus 16 --prefix ${savedir}/${flprefix}.king
fi

if [[ $1 == "divergence" ]]; then
	echo calculating divergence
	R CMD BATCH --vanilla "--args --prefix.in ${basedir}/beds/${flprefix}_pruned --file.seg ${savedir}/${flprefix}.king.seg --num_threads 16 --degree 4 --divThresh -0.02209709 --nRandomSNPs 0 --prefix.out ${savedir}/${flprefix}.divergence" ${basedir}/getDivergence_wrapper.R ${savedir}/getDivergence.Rout
fi

if [[ $1 == "unrelated" ]]; then
	echo extracting unrelated samples
	R CMD BATCH --vanilla "--args --prefix.in ${basedir}/beds/${flprefix}_pruned --file.seg ${savedir}/${flprefix}.king.seg --degree 4 --file.div ${savedir}/${flprefix}.divergence.div --prefix.out ${savedir}/${flprefix}_unrelated" ${basedir}/extractUnrelated_wrapper.R ${savedir}/extractUnrelated.Rout
fi

if [[ $1 == "pca" ]]; then
	echo running pca
	R  CMD BATCH --vanilla "--args --prefix.in ${basedir}/beds/${flprefix}_pruned --file.unrels ${savedir}/${flprefix}_unrelated.unrels --prefix.out ${savedir}/${flprefix}_pca --no_pcs 20 --num_threads 16 --no_iter 10" ${basedir}/runPCA_wrapper.R ${savedir}/runPCA.Rout
fi

if [[ $1 == "grm" ]]; then
	echo generating GRM
	R CMD BATCH --vanilla "--args --prefix.in ${basedir}/beds/${flprefix}_pruned --prefix.out ${savedir}/${flprefix}.sparseGRM --file.train ${savedir}/${flprefix}_unrelated.unrels --file.score ${savedir}/${flprefix}_pca.score --file.seg ${savedir}/${flprefix}.king.seg --num_threads 16 --no_pcs 20 --block.size 5000 --max.related.block 5000 --KINGformat.out FALSE --degree 4" ${basedir}/calcSparseGRM_wrapper.R ${savedir}/calcSparseGRM.Rout
fi

