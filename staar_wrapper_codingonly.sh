#!/bin/bash

# This script is submitted via the 'batch job' scripts
# and accepts the job id and which step of the pipeline to run as arguments
arrid=$2
cohort_analyze=discovery  # topmed  # discovery

# Path to this directory
scriptdir=/home/kd2630/STAARpipeline-MetaSTAAR-mod
source ${scriptdir}/makeenv2.sh
# Base directory where all output data will reside
basedir=/mnt/cuimc-med-pulm-general/STAAR/
#basedir=/home/kd2630/STAAR_10k/
# Specific sub-folder, pertaining to the current cohort being analyzed
#
if [[ $cohort_analyze == "discovery" ]]; then
	filesdir=/10k_Cohort_codingonly/
	subdir=""
	vcfdir=/mnt/cuimc-med-pulm-general/vcfs/10k_hg38/samples_to_keep/GENCODE_coding_subset_100bp_exon_flank/
	input_vcf_name_1=gencode.coding.curegn.wgs.freeze.2a.chr
	input_vcf_name_2=.hg38.filtered.samp.list.pruned.vcf.gz
	gds_file_name_1="10k_coding_chr_"
	gds_file_name_2=".gds"
elif [[ $cohort_analyze == "topmed" ]]; then
	filesdir=/TOPMed_codingonly/
	subdir=""
	vcfdir=/mnt/cuimc-med-pulm-general/vcfs/topmed_full_concat/GENCODE_coding_subset_100bp_exon_flank/
	input_vcf_name_1=gencode.coding.topmed_chr
	input_vcf_name_2=.vcf.gz
	gds_file_name_1="topmed_coding_chr_"
	gds_file_name_2=".gds"
else 
	echo "Error: incorrect cohort to analyze" >&2
	exit 1

fi
#filesdir=/10k_Cohort_rm_telo_qv/  #/TOPMed_Full_Cohort_grm/  #/10k_Cohort_rm_telo_qv/
# A separate output directory, Can be set to 'filesdir' if not running the MetaSTAAR analyses
savedir=/MetaSTAAR_Discovery_TOPMed_codingonly/ #/Gene_Noncoding/  #Indiv_Analysis/  #${filesdir}

# Where base STAAR files (gds, VariantInfo, etc.) are stored
dir_geno=${basedir}/${filesdir}
res_savedir=${basedir}/${savedir}${subdir}

# Parameters to some STAAR scripts:
# Allele frequency and count cutoffs,
# and variant type (set to 'variant' to include both SNVs and indels)
format_convert_to_gds=vcf
af_thresh=0.01
mac_cutoff=20
var_type=SNV  # SNV  # Indel

variant_output_path=${dir_geno}"/VariantInfo"
xsv="/home/kd2630/.cargo/bin/xsv"
# TODO this needs to change to the S3 location
#favor_path=${basedir}"/FAVORDB/"
# for testing with chr22
favor_path=/mnt/cuimc-med-pulm-general/tmp/

prestep_out=${dir_geno}/"/AssociationAnalysisPrestep/"

if [[ $cohort_analyze == "discovery" ]]; then
	# Including RV carriers
	phenofile=${dir_geno}/discovery_agesex_input_model_data.csv
	sgrmfile=${dir_geno}"/GRM/output/10k_chr_all.sparseGRM.sGRM.RData"
	null_model_fl=${dir_geno}/"staar_null_model/obj_nullmodel_All_w_qv.Rdata"
elif [[ $cohort_analyze == "topmed" ]]; then
	# Including RV carriers
	phenofile=${dir_geno}/topmed_agesex_input_model_data.csv
	sgrmfile=${dir_geno}"/GRM/output/topmed_chr_all.sparseGRM.sGRM.RData"
	null_model_fl=${dir_geno}/"staar_null_model/obj_nullmodel_w_qv.Rdata"
else 
	echo "Error: incorrect cohort to analyze" >&2
	exit 1
fi

# Names of output directories for each analysis
indivdir=/Individual_Variant_Analysis/
codingdir=/Gene_Centric_Coding_Analysis/
noncodingdir=/Gene_Centric_Noncoding_Analysis/
ncrnadir=/Gene_Centric_ncRNA_Analysis/
slidwindir=/Sliding_Window_Analysis/
dynamwindir=Dynamic_Window_Analysis/


echo For project $dir_geno #$savedir

# # Setup (These are run via 'run_staar_wrapper.sh')
# create gds files
# NOTE: this default script from STAARpipeline-Tutorial was created to convert vcf to gds and uses the seqVCF2GDS function
# from the SeqArray package. I have included the seqBED2GDS function from the same package to convert BED files from PLINK as well
if [[ $1 == "gds" ]]; then
	echo Creating GDS file for chrom $arrid
	inputvcf=${vcfdir}/${input_vcf_name_1}${arrid}${input_vcf_name_2}
	echo $inputvcf
	Rscript ${scriptdir}/STAARpipeline/convertVCF2GDS.R NULL ${format_convert_to_gds} ${dir_geno}/${gds_file_name_1}${arrid} 1 ${inputvcf} ${filesdir}
fi

# create variant list files
if [[ $1 == "listVariants" ]]; then
	echo Creating variant lists for chrom $arrid
	Rscript ${scriptdir}/STAARpipeline/FAVORannotator_csv/Varinfo_gds.R ${arrid} ${dir_geno} ${gds_file_name_1} ${gds_file_name_2} ${variant_output_path} ${scriptdir}/STAARpipeline/
fi

# annotate variants
if [[ $1 == "annotateVariants" ]]; then
	echo Annotating variants for chrom $arrid and including AlphaMissense
	Rscript ${scriptdir}/STAARpipeline/FAVORannotator_csv/Annotate_mod.R ${arrid} ${xsv} ${variant_output_path} ${favor_path} ${scriptdir}/STAARpipeline/
fi

# annotate gds files
if [[ $1 == "agds" ]]; then
	echo Creating aGDS files for chrom $arrid
	Rscript ${scriptdir}/STAARpipeline/FAVORannotator_csv/gds2agds.R ${arrid} ${dir_geno} ${gds_file_name_1} ${gds_file_name_2} ${variant_output_path}
fi

# run 'prestep' to create files with paths to input files
if [[ $1 == "prestep" ]]; then
	echo Running PreStep for
	Rscript ${scriptdir}/STAARpipeline/Association_Analysis_PreStep.r ${dir_geno} ${gds_file_name_1} ${gds_file_name_2} ${prestep_out}
fi

# create null model
if [[ $1 == "nullmodel" ]]; then
	echo Creating null model 
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Null_Model.r ${phenofile} ${null_model_fl} ${sgrmfile}
fi



# # Analyses (These are run via job submission using the corresponding script in 'STAARpipeline/batch jobs')
if [[ $1 == "individual_analysis" ]]; then
	echo Running individual variant analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Individual_Analysis.r ${arrid} ${res_savedir} ${indivdir} ${mac_cutoff} ${var_type} ${null_model_fl}
fi

if [[ $1 == "coding_analysis" ]]; then
	echo Running gene-centric coding analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_Coding.r ${arrid} ${dir_geno} ${res_savedir} ${codingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "noncoding_analysis" ]]; then
	echo Running gene-centric noncoding analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_Noncoding.r ${arrid} ${dir_geno} ${res_savedir} ${noncodingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "ncRNA_analysis" ]]; then
	echo Running gene-centric noncoding RNA analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_ncRNA.r ${arrid} ${dir_geno} ${res_savedir} ${ncrnadir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "sliding_window" ]]; then
	echo Running sliding window analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Sliding_Window.r ${arrid} ${res_savedir} ${slidwindir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "dynamic_window" ]]; then
	echo Running dynamic window analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Dynamic_Window.r ${arrid} ${res_savedir} ${dynamwindir} ${var_type} ${null_model_out}
fi


# # MetaSTAAR Analyses (These are run via job submission using the corresponding script in 'MetaSTAAR/batch jobs')
if [[ $1 == "individual_worker" ]]; then
	echo Running individual variant worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Individual_Analysis.r ${arrid} ${dir_geno} ${res_savedir} ${indivdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "coding_worker" ]]; then
	echo Running gene-centric coding worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_Coding.r ${arrid} ${dir_geno} ${res_savedir} ${codingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "coding_worker_null_mod_adj" ]]; then
	echo Running gene-centric coding worker analysis using gene-adjusted null models
	Rscript ${scriptdir}/MetaSTAAR/coding_MetaSTAARlite_worker_per_gene.R ${dir_geno} ${res_savedir} ${codingdir} ${af_thresh} ${var_type}
fi

if [[ $1 == "noncoding_worker" ]]; then
	echo Running gene-centric noncoding worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_Noncoding.r ${arrid} ${dir_geno} ${res_savedir} ${noncodingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "noncoding_worker_null_mod_adj" ]]; then
	echo Running gene-centric noncoding worker analysis using gene-adjusted null models
	Rscript ${scriptdir}/MetaSTAAR/noncoding_MetaSTAARlite_worker_per_gene.R ${dir_geno} ${res_savedir} ${noncodingdir} ${af_thresh} ${var_type}
fi

if [[ $1 == "ncRNA_worker" ]]; then
	echo Running gene-centric ncRNA worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_ncRNA.r ${arrid} ${dir_geno} ${res_savedir} ${ncrnadir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "individual_meta" ]]; then
	echo Running individual variant meta analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_Individual_Analysis.r ${arrid}
fi

if [[ $1 == "coding_meta" ]]; then
	echo Running gene-centric coding meta analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_Gene_Centric_Coding.r ${arrid}
fi

if [[ $1 == "noncoding_meta" ]]; then
	echo Running gene-centric noncoding meta analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_Gene_Centric_Noncoding.r ${arrid}
fi

if [[ $1 == "ncrna_meta" ]]; then
	echo Running gene-centric ncRNA meta analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_Gene_Centric_ncRNA.r ${arrid}
fi


# # MetaSTAAR cond analyses (These are run via job submission using the corresponding script in 'batch jobs')
if [[ $1 == "individual_worker_cond" ]]; then
	echo Running individual variant cond worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Individual_Analysis_cond.r ${dir_geno} ${res_savedir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "noncoding_worker_cond" ]]; then
	echo Running gene-centric noncoding cond worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_Noncoding_cond.r ${dir_geno} ${res_savedir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "ncRNA_worker_cond" ]]; then
	echo Running gene-centric ncRNA cond worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_ncRNA_cond.r ${dir_geno} ${res_savedir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "coding_worker_cond" ]]; then
	echo Running gene-centric coding cond worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_Coding_cond.r ${dir_geno} ${res_savedir} ${af_thresh} ${var_type} ${null_model_fl}
fi

