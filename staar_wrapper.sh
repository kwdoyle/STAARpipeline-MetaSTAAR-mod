# This script is submitted via the 'batch job' scripts
# and accepts the job id and which step of the pipeline to run as arguments
arrid=$2

# Path to this directory
scriptdir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/
source ${scriptdir}/makeenv2.sh
# Base directory where all output data will reside
basedir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/
# Specific sub-folder, pertaining to the current cohort being analyzed
filesdir=/10k_Cohort_rm_telo_qv/  #/TOPMed_Full_Cohort_grm/  #/10k_Cohort_rm_telo_qv/
# A separate output directory, Can be set to 'filesdir' if not running the MetaSTAAR analyses
savedir=/MetaSTAAR_Discovery_TOPMed/  #${filesdir}

# An optional subdirectory where to store/read outputs from.
subdir="/Discovery/"  #"/TOPMed/"  #"/Discovery/"  #""

# Where base STAAR files (gds, VariantInfo, etc.) are stored
dir_geno=${basedir}/${filesdir}
res_savedir=${basedir}/${savedir}${subdir} #${dir_geno}

# Parameters to some STAAR scripts:
# Allele frequency and count cutoffs,
# and variant type (set to 'variant' to include both SNVs and indels)
format_convert_to_gds=vcf
af_thresh=0.01
mac_cutoff=20
var_type=variant  # SNV  # Indel

gds_file_name_1="10k_chr_"  #"10k_chr_"  # "topmed_chr_"
gds_file_name_2=".gds"
variant_output_path=${dir_geno}"/VariantInfo"
xsv="/efs/garcia/users/kd2630/.cargo/bin/xsv"
favor_path=${basedir}"/FAVORDB/"

prestep_out=${dir_geno}/"/AssociationAnalysisPrestep/"

# Path to null model file (either to first create, or later use in analyses)
#null_model_fl=${dir_geno}/"staar_null_model/obj_nullmodel.Rdata"
# Use the 'all' model for discovery:
null_model_fl=${dir_geno}/"staar_null_model/obj_nullmodel_All_rm_qv.Rdata"
# For TOPMed:
# #null_model_fl=${dir_geno}/"staar_null_model/obj_nullmodel.Rdata"

# Names of output directories for each analysis
indivdir=/Individual_Variant_Analysis/
codingdir=/Gene_Centric_Coding_Analysis/
noncodingdir=/Gene_Centric_Noncoding_Analysis/
ncrnadir=/Gene_Centric_ncRNA_Analysis/
slidwindir=/Sliding_Window_Analysis/
dynamwindir=Dynamic_Window_Analysis/


echo For project $savedir

# # Setup (These are run via 'run_staar_wrapper.sh')
# create gds files
# NOTE: this default script from STAARpipeline-Tutorial was created to convert vcf to gds and uses the seqVCF2GDS function
# from the SeqArray package. I have included the seqBED2GDS function from the same package to convert BED files from PLINK as well
if [[ $1 == "gds" ]]; then
	echo Creating GDS file for chrom $arrid
	Rscript ${scriptdir}/STAARpipeline/convertVCF2GDS.R NULL ${format_convert_to_gds} ${dir_geno}/${savegdsprefix}${arrid} 1 ${vcfdir}/*chr${arrid}.*vcf.gz
fi

# create variant list files
if [[ $1 == "listVariants" ]]; then
	echo Creating variant lists for chrom $arrid
	Rscript ${scriptdir}/STAARpipeline/FAVORannotator_csv/Varinfo_gds.R ${arrid} ${dir_geno} ${gds_file_name_1} ${gds_file_name_2} ${variant_output_path} ${scriptdir}
fi

# annotate variants
if [[ $1 == "annotateVariants" ]]; then
	echo Annotating variants for chrom $arrid
	Rscript ${scriptdir}/STAARpipeline/FAVORannotator_csv/Annotate.R ${arrid} ${xsv} ${variant_output_path} ${favor_path} ${scriptdir}
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
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_Coding.r ${arrid} ${res_savedir} ${codingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "noncoding_analysis" ]]; then
	echo Running gene-centric noncoding analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_Noncoding.r ${arrid} ${res_savedir} ${noncodingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "ncRNA_analysis" ]]; then
	echo Running gene-centric noncoding RNA analysis for array id $arrid
	Rscript ${scriptdir}/STAARpipeline/STAARpipeline_Gene_Centric_ncRNA.r ${arrid} ${res_savedir} ${ncrnadir} ${af_thresh} ${var_type} ${null_model_fl}
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

if [[ $1 == "noncoding_worker" ]]; then
	echo Running gene-centric noncoding worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_Noncoding.r ${arrid} ${dir_geno} ${res_savedir} ${noncodingdir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "ncRNA_worker" ]]; then
	echo Running gene-centric ncRNA worker analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_worker_Gene_Centric_ncRNA.r ${arrid} ${dir_geno} ${res_savedir} ${ncrnadir} ${af_thresh} ${var_type} ${null_model_fl}
fi

if [[ $1 == "individual_meta" ]]; then
	echo Running individual variant meta analysis for array id $arrid
	Rscript ${scriptdir}/MetaSTAAR/MetaSTAARlite_Individual_Analysis.r ${arrid}
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

