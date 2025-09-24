#####################################################################
# Dynamic window analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(SCANG)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
arrayid <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[2]
resdir <- commandArgs(TRUE)[3]
#afthresh <- as.numeric(commandArgs(TRUE)[4])
variant_type <- commandArgs(TRUE)[4]
nullmodeldir <- commandArgs(TRUE)[5]

print(paste("Base dir is", basedir))
print(paste("Analysis dir is", resdir))
print(paste("Null model dir is", nullmodeldir))
print(paste("Using variant type:", variant_type))


## Number of jobs for each chromosome
jobs_num <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
#obj_nullmodel_SCANG_STAAR <- get(load(paste0(basedir, "/staar_null_model/obj_nullmodel_SCANG_STAAR.Rdata")))
obj_nullmodel_SCANG_STAAR <- get(load(paste0(nullmodeldir, "/obj_nullmodel_SCANG_STAAR.Rdata"))) 

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/Annotation_name_catalog.Rdata")))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
output_path <- paste0(basedir, "/", resdir)
print(paste("Will create directory:", output_path))
dir.create(output_path)

## output file name
#output_file_name <- "TOPMed_F5_LDL_SCANG"
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Dynamic_Window_Analysis")
print(paste("Output file name is", output_file_name))

###########################################################
#           Main Function 
###########################################################
## Number of jobs for SCANG
sum(jobs_num$scang_num)

chr <- which.max(arrayid <= cumsum(jobs_num$scang_num))
group.num <- jobs_num$scang_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$scang_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*1.5e6 + jobs_num$start_loc[chr]
end_loc <- min(start_loc + 1.5e6 - 1, jobs_num$end_loc[chr])

results_scang <- Dynamic_Window_SCANG(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel_SCANG_STAAR,
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

save(results_scang,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

