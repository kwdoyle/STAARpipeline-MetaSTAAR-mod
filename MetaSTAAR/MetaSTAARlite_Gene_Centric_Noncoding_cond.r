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
#file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond/",
#              "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond/")
#file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond_rs111521887/",
#              "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond_rs111521887/")
# Main:
#file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond_gwashits/",
#	      "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond_gwashits/")
# For DOCK2 conditioned on TERT gwas variants
#file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_TERT_gwas/",
#	      "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_TERT_gwas/")
# For DOCK2 conditioned on SPDL1 gwas variants
#file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_SPDL1_gwas/",
#	      "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_SPDL1_gwas/")
file.dir <- c("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Cond_codinghits/",
	      "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Cond_codinghits/")
file.prefix <- c("10k_Cohort_rm_telo_qv_Gene_Noncoding_Analysis_Cond", "TOPMed_Full_Cohort_grm_Gene_Noncoding_Analysis_Cond")
## Sample sizes of participating studies
sample.sizes <- c(3602, 6207)
## variant_type
variant_type <- "variant"
## cov_maf_cutoff
cov_maf_cutoff <- c(0.05, 0.05)

## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

## output path
#output_path <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis_Cond/"
#output_path <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis_Cond_gwashits/"
#output_path <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis_Cond_gwashits_DOCK2_on_TERT_gwas/"
#output_path <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis_Cond_gwashits_DOCK2_on_SPDL1_gwas/"
output_path <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis_Cond_codinghits/"
dir.create(output_path, recursive=T, showWarnings=F)
## output file name
output_file_name <- "noncodingMetaCond"

results_noncoding_cond <- c()

sumstat.file.path <- paste0(file.dir,file.prefix,"_sumstat.Rdata")
cov.file.path <- paste0(file.dir,file.prefix,"_cov.Rdata")
covcond.file.path <- paste0(file.dir,file.prefix,"_cov_cond.Rdata")
noncoding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)
noncoding_cov_cond_list <- sapply(covcond.file.path, function(x) mget(load(x)), simplify = TRUE)

# Perform the below for every significant hit
#sighits <- read.csv("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis/noncoding_sig.csv")
#sighits_filt <- sighits[sighits$MetaSTAAR.O < 3.57E-07 & sighits$X.SNV > 10, ]
#sighits_filt <- sighits_filt[order(sighits_filt$MetaSTAAR.O), ]
#sighits_nms <- unique(sighits_filt[, c("Gene.name", "Chr")])
## temporary mod to only analyze DOCK2
#sighits_nms <- sighits_nms[sighits_nms$Gene.name == "DOCK2", ]

sighits <- read.csv("/efs/garcia/users/kd2630/MetaSTAAR_noncoding_with_nearby_coding_hit.csv")
sighits_nms <- unique(sighits[ ,c("Gene.name", "Chr")])

#chr <- 11
#gene_name <- "TOLLIP"

for (i in 1:nrow(sighits_nms)) {
	chr <- sighits_nms$Chr[i] #11 #19
	gene_name <- sighits_nms$Gene.name[i]
	print(paste("Chr", chr))
	print(paste("Gene", gene_name))

	# NOTE: the appended names for these got changed to match the variable names from when I manually combined all the genes from the worker script into one list
	noncoding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
	  noncoding_sumstat_list[[paste0(x,".sumstat_dat")]][[gene_name]]  # sumstat_dat  # noncoding_sumstat
	})
	noncoding_cov_gene_list <- lapply(cov.file.path, function(x) {
	  noncoding_cov_list[[paste0(x,".cov_dat")]][[gene_name]] # cov_dat  # noncoding_cov
	})
	noncoding_cov_cond_gene_list <- lapply(covcond.file.path, function(x) {
	  noncoding_cov_cond_list[[paste0(x,".covcond_dat")]][[gene_name]]  # covcond_dat  # noncoding_cov_cond
	})

	if ( any(unlist(lapply(noncoding_cov_cond_gene_list, length)) == 0) ) {
		message("Gene has no cov cond output")
		next
		
	}
	results_cond <- noncoding_MetaSTAARlite_cond(chr=chr,gene_name=gene_name,
						     sample.sizes=sample.sizes,noncoding_sumstat_gene_list=noncoding_sumstat_gene_list,
						     noncoding_cov_gene_list=noncoding_cov_gene_list,
						     noncoding_cov_cond_gene_list=noncoding_cov_cond_gene_list,
						     cov_maf_cutoff=cov_maf_cutoff,
						     rare_maf_cutoff=0.01,rv_num_cutoff=2,
						     check_qc_label=TRUE,variant_type=variant_type,
						     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	# put into list with name of current gene
	results_cond_l <- list(results_cond)
	names(results_cond_l) <- gene_name

	results_noncoding_cond <- append(results_noncoding_cond, results_cond_l)

}

#save(results_noncoding_cond, file=paste0(output_path, output_file_name, "_", gene_name, ".Rdata"))
save(results_noncoding_cond, file=paste0(output_path, output_file_name, ".Rdata"))

