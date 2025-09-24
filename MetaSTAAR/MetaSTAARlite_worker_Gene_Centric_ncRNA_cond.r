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

## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"

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

known_loci <- read.table("~/noncoding_telo/STAAR/gwas_hits_shared_disc_tm.tsv")
names(known_loci) <- c("CHR", "POS", "REF", "ALT")

# Perform the below for every significant hit
sighits <- read.csv("~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Noncoding_Meta_Analysis/noncoding_sig.csv")
sighits_filt <- sighits[sighits$MetaSTAAR.O < 3.57E-07 & sighits$X.SNV > 10 & sighits$Category == "ncRNA", ]
sighits_filt <- sighits_filt[order(sighits_filt$MetaSTAAR.O), ]
sighits_nms <- unique(sighits_filt[, c("Gene.name", "Chr")])

## output path
#output_path <- "/path_to_the_output_file/"
### output file name
#output_file_name <- "JHS_MIR4497_ncRNA"
outdirnm <- "/Gene_Centric_ncRNA_Cond_gwashits/" 

output_path <- paste0(resdir, "/", outdirnm)
print(paste("Will create directory:", output_path))
dir.create(output_path, recursive=TRUE)
## output file name
pthsplt <- unlist(strsplit(basedir, "\\/"))
nmuse <- pthsplt[length(pthsplt)]
output_file_name <- paste0(nmuse, "_Gene_ncRNA_Analysis_Cond")

print(paste("Output file name is", output_file_name))


###############################
# ncRNA
###############################
# Loop over all sig, ncRNA hits
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

	ncRNA_sumstat <- list()
	ncRNA_cov <- list()
	ncRNA_cov_cond <- list()

	#known_loci_MIR4497 <- data.frame(CHR=12,POS=109543379,REF="C",ALT="CG")
	#gene_name <- "MIR4497"
	results_temp <- ncRNA_MetaSTAARlite_worker(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
						   known_loci=known_loci_chr,
						   cov_maf_cutoff=0.05,signif.digits=NULL,
						   QC_label=QC_label,check_qc_label=TRUE,variant_type=variant_type,
						   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
						   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	ncRNA_sumstat[[gene_name]] <- results_temp$summary_stat
	ncRNA_cov[[gene_name]] <- results_temp$GTSinvG_rare
	ncRNA_cov_cond[[gene_name]] <- results_temp$cov_cond

	save(ncRNA_sumstat, file=paste0(output_path, output_file_name, "_", gene_name, "_sumstat.Rdata"), compress = "xz")
	save(ncRNA_cov, file=paste0(output_path, output_file_name, "_", gene_name, "_cov.Rdata"), compress = "xz")
	save(ncRNA_cov_cond, file=paste0(output_path, output_file_name, "_", gene_name, "_cov_cond.Rdata"), compress = "xz")

	seqClose(genofile)
}

