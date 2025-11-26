#!/bin/bash

scriptdir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/
basedir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/
# directory with pheno files
#indir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models
indir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models_w_c1q
#indir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/coding_gene_vnt_null_models
#indir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/coding_gene_vnt_null_models_w_c1q

#sgrmfile=${dir_geno}"/GRM/output/10k_chr_all.sparseGRM.sGRM.RData"
#sgrmfile=${dir_geno}"/GRM/output/topmed_chr_all.sparseGRM.sGRM.RData"

# loop through them and make corresponding null model
for file in "$indir"/*.csv; do
	# base name
	#fname=$(basename "$file")
	phenofile=$file
	echo Phenotype File:
	echo $phenofile
	#file $phenofile
	echo 

	# use appropriate sgrm file
	if [[ "$file" == *"discovery"* ]]; then
		sgrmfile=${basedir}/"10k_Cohort_rm_telo_qv/GRM/output/10k_chr_all.sparseGRM.sGRM.RData"
	elif [[ "$file" == *"topmed"* ]]; then
		sgrmfile=${basedir}"/TOPMed_Full_Cohort_grm/GRM/output/topmed_chr_all.sparseGRM.sGRM.RData"
	fi

	echo GRM File:
	echo $sgrmfile
	#file $sgrmfile
	echo
	
	# make null model file name
	null_model_fl="${phenofile%.csv}_nullmodel.Rdata"
	echo Null Model File to Save:
	echo $null_model_fl
	echo

	# run null model staar script
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Null_Model.r ${phenofile} ${null_model_fl} ${sgrmfile}

done

