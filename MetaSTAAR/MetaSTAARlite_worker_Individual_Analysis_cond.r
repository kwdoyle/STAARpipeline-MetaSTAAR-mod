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
library(dplyr)

###########################################################
#           User Input
###########################################################
#arrayid <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[1]
resdir <- commandArgs(TRUE)[2]
#outdirnm <- commandArgs(TRUE)[4]
afthresh <- as.numeric(commandArgs(TRUE)[3])
variant_type <- commandArgs(TRUE)[4]
nullmodelfl <- commandArgs(TRUE)[5]

print(paste("Base dir is", basedir))
print(paste("Analysis dir is", resdir))
print(paste("Null model file is", nullmodelfl))
print(paste("Using MAF cutoff of", afthresh))
print(paste("Using variant type:", variant_type))

alpha <- 5E-08

## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/Annotation_name_catalog.Rdata")))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- FALSE
## Annotation name
Annotation_name <- NULL

known_loci <- read.table("~/noncoding_telo/STAAR/gwas_hits_shared_disc_tm.tsv")
names(known_loci) <- c("CHR", "POS", "REF", "ALT")

sighits <- read.csv("~/noncoding_telo/STAAR/variants_to_cond.csv")
#sighits_filt <- sighits[sighits$pvalue < alpha, ]
sighits_range <- sighits %>%
	group_by(CHR, Gene) %>%
	summarise(from = min(POS, na.rm=T),
		  to = max(POS, na.rm=T))

# temporary mod to only analyze DOCK2
gwas_gene_snps_use <- "SPDL1"  #"SPDL1"  #"TERT"
sighits_range <- sighits_range[sighits_range$Gene == "DOCK2", ]

#outdirnm <- "/Individual_Variant_Analysis_Cond_gwashits/" 
#outdirnm <- "/Individual_Variant_Analysis_Cond_gwashits_DOCK2_on_TERT_gwas//" 
outdirnm <- "/Individual_Variant_Analysis_Cond_gwashits_DOCK2_on_SPDL1_gwas//" 

output_path <- paste0(resdir, "/", outdirnm)
print(paste("Will create directory:", output_path))
dir.create(output_path, recursive=TRUE)
## output file name
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Individual_Variant_Analysis_Cond")
print(paste("Output file name is", output_file_name))

###############################
# LDLR individual analysis
###############################
for (i in 1:nrow(sighits_range)) {

	chr <- sighits_range$CHR[i] #11 #19
	gene_name <- sighits_range$Gene[i]
	print(paste("Chr", chr))
	print(paste("Gene", gene_name))

	## aGDS file
	agds.path <- agds_dir[chr]
	genofile <- seqOpen(agds.path)

	# Subset known_loci for the current chr
	known_loci_chr <- known_loci[known_loci$CHR == chr, ]
	# temp mod to condition only on TERT and SPDL1 gwas hits separately
	if (gwas_gene_snps_use == "TERT") {
		known_loci_chr <- known_loci_chr[dplyr::between(known_loci_chr$POS, 1200000, 1300000), ]
	} else if (gwas_gene_snps_use == "SPDL1") { 
		known_loci_chr <- known_loci_chr[known_loci_chr$POS > 100000000, ]
	}

	print("Conditioning on variants:")
	print(known_loci_chr)

	start_loc <- sighits_range$from[i]
	end_loc <- sighits_range$to[i]
	results_temp <- individual_analysis_MetaSTAARlite_worker(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel,
								 known_loci=known_loci_chr,
								 subsegment.size=5e4,
								 QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
								 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
								 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	individual_analysis_sumstat <- results_temp$summary_stat
	individual_analysis_cov_cond <- results_temp$cov_cond

	save(individual_analysis_sumstat, file=paste0(output_path, output_file_name, "_", gene_name, "_sumstat.Rdata"), compress = "xz")
	save(individual_analysis_cov_cond, file=paste0(output_path, output_file_name, "_", gene_name,"_cov_cond.Rdata"), compress = "xz")

	seqClose(genofile)
}

