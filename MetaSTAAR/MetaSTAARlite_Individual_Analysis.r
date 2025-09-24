rm(list=ls())
gc()

## load required packages
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## Directories of the study-specific summary statistics file folders
file.dir <- c("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Individual_Variant_Analysis/",
              "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Individual_Variant_Analysis/")
file.prefix <- c("10k_Cohort_rm_telo_qv_Individual_Variant_Analysis", "TOPMed_Full_Cohort_grm_Individual_Variant_Analysis")
## Sample sizes of participating studies
sample.sizes <- c(3602, 6207)

## output path
output_path <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Individual_Meta_Analysis/"
dir.create(output_path, recursive=T, showWarnings=F)
## output file name
output_file_name <- "individualMeta"
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat_",arrayid,".Rdata")
individual_analysis_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)

results_individual_analysis <- individual_analysis_MetaSTAARlite(sample.sizes=sample.sizes,
                                                                 sumstat.list=individual_analysis_sumstat_list,
                                                                 mac_cutoff=20,check_qc_label=TRUE)

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

