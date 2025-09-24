#####################################################################
# Sliding window analysis using STAARpipeline
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
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## input array id from batch file (Harvard FAS RC cluster)
# note: actually use the intended large number of array jobs,
# splitting each chromosome up into many jobs.
# this analysis will likely take too long, otherwise
arrayid <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[2]
resdir <- commandArgs(TRUE)[3]
afthresh <- as.numeric(commandArgs(TRUE)[4])
variant_type <- commandArgs(TRUE)[5]
#nullmodeldir <- commandArgs(TRUE)[6]
nullmodelfl <- commandArgs(TRUE)[6]

print(paste("Base dir is", basedir))
print(paste("Analysis dir is", resdir))
print(paste("Null model file is", nullmodelfl))
print(paste("Using MAF cutoff of", afthresh))
print(paste("Using variant type:", variant_type))

## Number of jobs for each chromosome
jobs_num <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
#obj_nullmodel <- get(load(paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata")))
#obj_nullmodel <- get(load(paste0(nullmodeldir, "/obj_nullmodel.Rdata"))) 
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## sliding_window_length
sliding_window_length <- 2000

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
#output_path <- paste0(basedir, "/Sliding_Window_Analysis/")
output_path <- paste0(basedir, "/", resdir)
print(paste("Will create directory:", output_path))
dir.create(output_path)

## output file name
#output_file_name <- "TOPMed_F5_LDL_Sliding_Window"
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Sliding_Window_Analysis")
print(paste("Output file name is", output_file_name))

###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$sliding_window_num))
group.num <- jobs_num$sliding_window_num[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(jobs_num$sliding_window_num)[chr-1]
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

start_loc <- (groupid-1)*5e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + (sliding_window_length/2)*20 - 1

results_sliding_window <- c()
for(kk in 1:(5e6/((sliding_window_length/2)*20)))
{
  print(kk)
  start_loc_sub <- start_loc + (sliding_window_length/2)*20*(kk-1)
  end_loc_sub <- end_loc + (sliding_window_length/2)*20*(kk-1) + (sliding_window_length/2)
  
  end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[chr])
  
  results <- c()
  if(start_loc_sub < end_loc_sub)
  {
    results <- try(Sliding_Window(chr=chr,start_loc=start_loc_sub,end_loc=end_loc_sub,
                                  sliding_window_length=sliding_window_length,type="multiple",
                                  genofile=genofile,obj_nullmodel=obj_nullmodel,
                                  rare_maf_cutoff=afthresh,rv_num_cutoff=2,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    
    if(class(results)[1]!="try-error")
    {
      results_sliding_window <- rbind(results_sliding_window,results)
    }
    
  }
}

save(results_sliding_window,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)

