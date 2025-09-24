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

# set to FALSE if want variant statistics from before conditioning
# NOTE: setting this to FALSE actually does nothing for me, as the sumstat output is the same regardless if known_loci are provided.
# the only difference is whether cov_cond is generated.
calculate_conditional <- TRUE 

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
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

### Correct over all gwas hits, per gene
#gwas_hits <- read.csv("~/noncoding_telo/STAAR/Published GWAS hits - lifted - posrefalt filledin.csv")
## filter for those only with ref and alt
#known_loci <- unique(gwas_hits[!is.na(gwas_hits$true_major), c("CHR", "POS", "true_major", "true_minor")])
#names(known_loci)[3:4] <- c("REF", "ALT")
# Use only gwas hits present in both cohorts
known_loci <- read.table("~/noncoding_telo/STAAR/gwas_hits_shared_disc_tm.tsv")
names(known_loci) <- c("CHR", "POS", "REF", "ALT")

#known_loci_11 <- known_loci[known_loci$CHR == 11, ]
#known_loci_9 <- known_loci[known_loci$CHR == 9, ]

# Perform the below for every significant hit
sighits <- read.csv("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis/noncoding_sig.csv")
sighits_filt <- sighits[sighits$MetaSTAAR.O < 3.57E-07 & sighits$X.SNV > 10, ]
sighits_filt <- sighits_filt[order(sighits_filt$MetaSTAAR.O), ]
sighits_nms <- unique(sighits_filt[, c("Gene.name", "Chr")])
# temporary mod to only analyze DOCK2
gwas_gene_snps_use <- "TERT"  #"SPDL1"  #"TERT"
sighits_nms <- sighits_nms[sighits_nms$Gene.name == "DOCK2", ]

### output path
#output_path <- "/path_to_the_output_file/"
### output file name
#output_file_name <- "JHS_LDLR_noncoding"

## output path
#outdirnm <- "/Gene_Centric_Noncoding_Cond/"  # 1219991 rs35705950
#outdirnm <- "/Gene_Centric_Noncoding_Cond_rs111521887/" # 1291476
#outdirnm <- "/Gene_Centric_Noncoding_Cond_gwashits/" 
#outdirnm <- "/Gene_Centric_Noncoding_NotCond_gwashits/" 
# For DOCK2 separately for TERT and SPDL1 variants
outdirnm <- "/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_TERT_gwas/" 
#outdirnm <- "/Gene_Centric_Noncoding_Cond_gwashits_DOCK2_on_SPDL1_gwas/" 

output_path <- paste0(resdir, "/", outdirnm)
print(paste("Will create directory:", output_path))
dir.create(output_path, recursive=TRUE)
## output file name
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
if (calculate_conditional) {
	output_file_name <- paste0(nmuse, "_Gene_Noncoding_Analysis_Cond")
} else {
	output_file_name <- paste0(nmuse, "_Gene_Noncoding_Analysis_NotCond")
}

print(paste("Output file name is", output_file_name))

###############################
# LDLR noncoding
###############################
# Loop over all sig hits
for (i in 1:nrow(sighits_nms)) {

	chr <- sighits_nms$Chr[i] #11 #19
	gene_name <- sighits_nms$Gene.name[i]
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

	noncoding_sumstat <- list()
	noncoding_cov <- list()
	noncoding_cov_cond <- list()

	if (!calculate_conditional) {
		print("Calculating unconditional output; setting known_loci to NULL")
		known_loci_chr <- NULL
	}

	# For now, just correct for each snp at a time
	#known_loci <- data.frame(CHR=11, POS=1219991, REF="G", ALT="T")
	#known_loci <- data.frame(CHR=11, POS=1291476, REF="C", ALT="G")

	#gene_name <- "TOLLIP"  #"LDLR"


	results_temp <- noncoding_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
						       known_loci=known_loci_chr,
						       cov_maf_cutoff=0.05,signif.digits=NULL,
						       QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
						       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	noncoding_sumstat[[gene_name]] <- results_temp$summary_stat_list
	noncoding_cov[[gene_name]] <- results_temp$GTSinvG_rare_list
	noncoding_cov_cond[[gene_name]] <- results_temp$cov_cond_list

	save(noncoding_sumstat, file=paste0(output_path, output_file_name, "_", gene_name, "_sumstat.Rdata"), compress = "xz")
	save(noncoding_cov, file=paste0(output_path, output_file_name, "_", gene_name, "_cov.Rdata"), compress = "xz")
	save(noncoding_cov_cond, file=paste0(output_path, output_file_name, "_", gene_name, "_cov_cond.Rdata"), compress = "xz")

	seqClose(genofile)
}

