library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

# TODO set the params from main bash script and run and save for discovery and topmed
basedir <- commandArgs(TRUE)[1]
resdir <- commandArgs(TRUE)[2]
outdirnm <- commandArgs(TRUE)[3]
afthresh <- as.numeric(commandArgs(TRUE)[4])
variant_type <- commandArgs(TRUE)[5]
#nullmodeldir <- commandArgs(TRUE)[6]
# just hard code the null model dir in these scripts.
#nullmodeldir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models/"
# for conditioning C1QB synonymous on C1QA and C1QC 
nullmodeldir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models_w_c1q/"

# For testing: (eventually will be passed to this script like the others in staar_wrapper.sh)
#basedir="~/noncoding_telo/STAAR//10k_Cohort_rm_telo_qv/"
#resdir="~/noncoding_telo/STAAR//MetaSTAAR_Discovery_TOPMed//Discovery/"
#afthresh=0.01
#variant_type="variant"
## TODO this will eventually instead point to the directory with all the null models pertaining to whichever
## model is adjusted for the current gene of interest
## (e.g., if adjusting the coding hit in GOT2, then this will be the null model with the params for 'has_variant_in_GOT2_[noncoding category]'
##nullmodelfl="~/noncoding_telo/STAAR//10k_Cohort_rm_telo_qv/staar_null_model/obj_nullmodel_All_rm_qv.Rdata"
#nullmodeldir="/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models/"

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

# eventually use all genes from here with in_coding_and_noncoding = TRUE
# note: not sure how sheet 2 is different; don't know what christine did to make this sheet
#hits <- readxl::read_xlsx("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/MetaSTAAR_discovery_topmed_noncoding_and_coding_all_hits_metastaar_o_less_0_05_ckg.xlsx", sheet=1)
#hits <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/16_coding_hits_w_noncoding.csv")
hits <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/18_coding_hits_w_noncoding_w_c1q.csv")
# this has 795 genes. Will need to limit this down somehow (whatever was done to make sheet 2?)
#hits_use <- hits[hits$in_coding_and_noncoding == TRUE, ]
# for now, just select the main ones I'm working with
#genes_chk_tmp <- c("GOT2", "SAMD9L", "ARL2", "ARL2-SNX15", "P3H3")
#hits_use <- hits_use[hits_use$`Gene name` %in% c("GOT2", "SAMD9L", "ARL2", "ARL2-SNX15", "P3H3"), ]
#hits_use <- hits_use[grep(paste(genes_chk_tmp, collapse="|"), hits_use$`Gene name`), ]
hits_use <- hits

# only need gene and chr
# Note: ALR2-SNX15 isn't tested here b/c specifically for SNX15, there is no corresponding coding and noncoding hit together.
hits_use_tbl <- unique(hits_use[,c("Gene.name", "Chr")])
# add a second C1QB so that I cond on C1QA and C1QC
hits_use_tbl <- rbind(hits_use_tbl, hits_use_tbl[which(hits_use_tbl$Gene.name == "C1QB"), ])


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
output_file_name <- paste0(nmuse, "_Gene_Coding_Analysis")
print(paste("Output file name is", output_file_name))


# # Main function
coding_sumstat <- list()
coding_cov <- list()
# indicator for if I analyzed the first C1QB
ran_1st_c1qb <- FALSE
for (i in 1:nrow(hits_use_tbl)) {
	chr <- hits_use_tbl$Chr[i]
	gene_name <- hits_use_tbl$Gene.name[i]
	print(paste(gene_name, "Chr:",  chr, sep=" "))
	if (gene_name == "C1QB") {
		if (ran_1st_c1qb) {
			print("Cond on C1QC")
			nm_gn_use <- "C1QC"
			savenm <- paste(gene_name, "on", nm_gn_use, sep="_")
		} else {
			print("Cond on C1QA")
			nm_gn_use <- "C1QA"
			savenm <- paste(gene_name, "on", nm_gn_use, sep="_")
			ran_1st_c1qb <- TRUE
		}
	} else if (gene_name %in% c("C1QA", "C1QC")) {
	   	print("Skipping")
		next
	} else {
		nm_gn_use <- gene_name
		savenm <- gene_name
	}

#	} else if (gene_name == "P3H3") {
#		print("Cond on CDCA3")
#		nm_gn_use <- "CDCA3"
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

	results_temp <- coding_MetaSTAARlite_worker(chr=chr, gene_name=gene_name, genofile=genofile, obj_nullmodel=obj_nullmodel, genes=genes_info,
					      cov_maf_cutoff=0.05, signif.digits=NULL,
					      QC_label=QC_label, check_qc_label=TRUE, variant_type=variant_type,
					      Annotation_dir=Annotation_dir, Annotation_name_catalog=Annotation_name_catalog,
					      Use_annotation_weights=Use_annotation_weights, Annotation_name=Annotation_name)

	coding_sumstat[[savenm]] <- results_temp$summary_stat_list
	coding_cov[[savenm]] <- results_temp$GTSinvG_rare_list

	seqClose(genofile)
}

# I think I can save the entire output accross genes as a single file 
save(coding_sumstat, file=paste0(output_path, output_file_name, "_sumstat", ".Rdata"), compress = "xz")
save(coding_cov, file=paste0(output_path, output_file_name, "_cov", ".Rdata"), compress = "xz")

