###########################################################
# Individual analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## input array id from batch file (Harvard FAS RC cluster)
arrayid <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[2]
resdir <- commandArgs(TRUE)[3]
mac_cutoff <- as.numeric(commandArgs(TRUE)[4])
variant_type <- commandArgs(TRUE)[5]
nullmodelfl <- commandArgs(TRUE)[6]

print(paste("Base dir is", basedir))
print(paste("Analysis dir is", resdir))
print(paste("Null model file is", nullmodelfl))
print(paste("Using MAC cutoff of", mac_cutoff))
print(paste("Using variant type:", variant_type))

## Number of jobs for each chromosome
jobs_num <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

## QC_label
QC_label <- "annotation/filter"
## variant_type
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
#output_path <- paste0(basedir, "/Individual_Variant_Analysis/")
output_path <- paste0(basedir, "/", resdir)
print(paste("Will create directory:", output_path))
dir.create(output_path)
## output file name
# # use last part of basedir in output file name
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Individual_Analysis")
print(paste("Output file name is", output_file_name))

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
group.num <- jobs_num$individual_analysis_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 10e6 - 1
end_loc <- min(end_loc,jobs_num$end_loc[chr])

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=mac_cutoff, # 20,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

