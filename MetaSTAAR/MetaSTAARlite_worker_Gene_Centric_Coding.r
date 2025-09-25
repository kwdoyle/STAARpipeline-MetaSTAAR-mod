rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## input array id from batch file
arrayid <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[2]
resdir <- commandArgs(TRUE)[3]
outdirnm <- commandArgs(TRUE)[4]
afthresh <- as.numeric(commandArgs(TRUE)[5])
variant_type <- commandArgs(TRUE)[6]
nullmodelfl <- commandArgs(TRUE)[7]

print(paste("Base dir is", basedir))
print(paste("Analysis dir is", resdir))
print(paste("Null model file is", nullmodelfl))
print(paste("Using MAF cutoff of", afthresh))
print(paste("Using variant type:", variant_type))

## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

## QC_label
QC_label <- "annotation/filter"

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
output_path <- paste0(resdir, "/", outdirnm)
print(paste("Will create directory:", output_path))
dir.create(output_path, recursive=TRUE)
## output file name
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Gene_Coding_Analysis")
print(paste("Output file name is", output_file_name))

###########################################################
#           Main Function 
###########################################################
## gene number in job
gene_num_in_array <- 50
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))
group.num <- group.num.allchr[chr]

if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}

## exclude large coding masks
if(arrayid==57)
{
  sub_seq_id <- setdiff(sub_seq_id,840)
}

if(arrayid==112)
{
  sub_seq_id <- setdiff(sub_seq_id,c(543,544))
}

if(arrayid==113)
{
  sub_seq_id <- setdiff(sub_seq_id,c(575,576,577,578,579,580,582))
}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

genes <- genes_info

coding_sumstat <- list()
coding_cov <- list()
for(kk in sub_seq_id)
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results_temp <- coding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,genes=genes,
                                              cov_maf_cutoff=0.05,signif.digits=NULL,
                                              QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  coding_sumstat[[gene_name]] <- results_temp$summary_stat_list
  coding_cov[[gene_name]] <- results_temp$GTSinvG_rare_list
}

save(coding_sumstat,file=paste0(output_path,output_file_name,"_sumstat_",arrayid,".Rdata"),compress = "xz")
save(coding_cov,file=paste0(output_path,output_file_name,"_cov_",arrayid,".Rdata"),compress = "xz")

seqClose(genofile)

