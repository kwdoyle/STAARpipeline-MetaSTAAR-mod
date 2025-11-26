library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

basedir <- commandArgs(TRUE)[1]
resdir <- commandArgs(TRUE)[2]
outdirnm <- commandArgs(TRUE)[3]
afthresh <- as.numeric(commandArgs(TRUE)[4])
variant_type <- commandArgs(TRUE)[5]
#nullmodeldir <- commandArgs(TRUE)[6]
# just hard code the null model dir in these scripts.
#nullmodeldir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/coding_gene_vnt_null_models/"
# for conditioning C1QA and C1QC on C1QB synonymous
nullmodeldir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/coding_gene_vnt_null_models_w_c1q/"


nullmodelfls <- list.files(nullmodeldir, pattern="*nullmodel.Rdata", full.names=T)
# subset for current cohort 
if (grepl("discovery", basename(resdir), ignore.case=T)) {
	print("USING DISCOVERY")
	nullmodelfls <- nullmodelfls[grep("discovery", basename(nullmodelfls), ignore.case=T)]
} else if (grepl("topmed", basename(resdir), ignore.case=T)) {
	print("USING TOPMED")
	nullmodelfls <- nullmodelfls[grep("topmed", basename(nullmodelfls), ignore.case=T)]
} else {
	stop("Cohort could not be determined in order to select null models to use.")
}


#hits <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/16_coding_hits_w_noncoding.csv")
# note: when C1QA and C1QB are analyzed, use the C1QB coding has_vnt data
# and skip C1QB as a gene to analyze by itself
hits <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/18_coding_hits_w_noncoding_w_c1q.csv")
hits_use <- hits
# only need gene and chr
# Note: ALR2-SNX15 isn't tested here b/c specifically for SNX15, there is no corresponding coding and noncoding hit together.
hits_use_tbl <- unique(hits_use[,c("Gene.name", "Chr")])


agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
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
output_file_name <- paste0(nmuse, "_Gene_Noncoding_Analysis")
print(paste("Output file name is", output_file_name))


# # Main function
noncoding_sumstat <- list()
noncoding_cov <- list()
for (i in 1:nrow(hits_use_tbl)) {
	chr <- hits_use_tbl$Chr[i]
	gene_name <- hits_use_tbl$Gene.name[i]
	print(paste(gene_name, "Chr:",  chr, sep=" "))
	if (gene_name == "C1QB") {
		print("Skipping.")
		next
	} else if (gene_name %in% c("C1QA", "C1QC")) {
	   	print("Cond on C1QB")
		nm_gn_use <- "C1QB"
	} else if (gene_name == "P3H3") {
		print("Cond on CDCA3")
		nm_gn_use <- "CDCA3"
	} else {
		nm_gn_use <- gene_name
	}

	# genes_info is a built-in data frame from STAAR
	#genes_info_chr <- genes_info[genes_info[,2]==chr,]
	# can just subset for current gene instead
	genes_info_gn <- genes_info[genes_info$hgnc_symbol == gene_name, ]

	# get null model for current gene
	null_model_gn <- grep(paste0("_", nm_gn_use, "_"), nullmodelfls, value=T)
	print(paste("Using null model:", null_model_gn))
	obj_nullmodel <- get(load(paste0(null_model_gn))) 

	## aGDS file
	agds.path <- agds_dir[chr]
	genofile <- seqOpen(agds.path)

	results_temp <- noncoding_MetaSTAARlite_worker(chr=chr, gene_name=gene_name, genofile=genofile, obj_nullmodel=obj_nullmodel, 
					      cov_maf_cutoff=0.05, signif.digits=NULL,
					      QC_label=QC_label, check_qc_label=TRUE, variant_type=variant_type,
					      Annotation_dir=Annotation_dir, Annotation_name_catalog=Annotation_name_catalog,
					      Use_annotation_weights=Use_annotation_weights, Annotation_name=Annotation_name)

	noncoding_sumstat[[gene_name]] <- results_temp$summary_stat_list
	noncoding_cov[[gene_name]] <- results_temp$GTSinvG_rare_list

	seqClose(genofile)
}

# I think I can save the entire output accross genes as a single file 
save(noncoding_sumstat, file=paste0(output_path, output_file_name, "_sumstat", ".Rdata"), compress = "xz")
save(noncoding_cov, file=paste0(output_path, output_file_name, "_cov", ".Rdata"), compress = "xz")

