library(dplyr)

# TODO the has_noncoding_variants_[cohort].csv files come from gt_variant_presence_absense.R
# and will eventually need to be updated to use GT data from ALL genes
# from MetaSTAAR_discovery_topmed_noncoding_and_coding_all_hits_metastaar_o_less_0_05_ckg.xlsx
# where in_coding_and_noncoding is TRUE

to_process <- "coding"  # "noncoding"
cohort <- "topmed"  # "discovery" "topmed"

if (cohort == "discovery") {
	message("processing discovery")
	phenfl <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR//10k_Cohort_rm_telo_qv/discovery_cohort_844_controls_rm_telo_gene_carriers_input_model_data.csv"
	vntcols <- read.csv(paste0("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_", to_process, "_variants_discovery.csv"), row.names=1)

} else if (cohort == "topmed") {
	message("processing topmed")
	phenfl <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR//TOPMed_Full_Cohort_grm/topmed_input_model_data_rv_rm.csv"
	vntcols <- read.csv(paste0("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_", to_process, "_variants_topmed.csv"), row.names=1)
} else {
	stop("incorrect cohort specified.")
}

savedir <- paste0("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/", to_process, "_gene_vnt_null_models/")

#vntcols <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_noncoding_variants_discovery.csv", row.names=1)
#vntcols <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_coding_variants_discovery.csv", row.names=1)

#vntcols <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_noncoding_variants_topmed.csv", row.names=1)
#vntcols <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/has_coding_variants_topmed.csv", row.names=1)

#savedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/noncoding_gene_vnt_null_models/"
#savedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/coding_gene_vnt_null_models/"
dir.create(savedir, showWarnings=F)

phen <- read.csv(phenfl)

#genes_include <- c("GOT2", "SAMD9L", "C1QB", "ARL2", "ARL2-SNX15", "RCC1L", "CDCA3", "P3H3")
# these are all the noncoding hits with a corresponding coding hit.
#genes_include <- c("GOT2", "SAMD9L", "ARL2", "ARL2-SNX15", "P3H3")
# get all genes to include from the vntcol names
colnms <- names(vntcols)
genes_include <- sub("^has_vnt_in_([A-Za-z0-9]+).*", "\\1", colnms)

# save new pheno files per gene
for (gn in genes_include) {
	print(gn)
	gn_split <- unlist(strsplit(gn, "-"))
	#vntcols_use <- vntcols[, grep(gn, names(vntcols)), drop=F]
	#vntcols_use <- vntcols[, grep(paste(genes_include, collapse="|"), names(vntcols))]

	# match names the other way around so that ARL2-SNX15 will include the ARL2 variant columns
	cols_matched <- grep(paste(gn_split, collapse="|"), colnms, value=TRUE)
#	matched_cols <- colnms[
#			  sapply(colnms, function(col) any(sapply(gn, function(g) grepl(g, col))))
#			]
	vntcols_use <- vntcols[, cols_matched, drop=F]

	vntcols_use$sample.id <- rownames(vntcols_use)

	phen_w_vntcols <- phen %>% left_join(vntcols_use, by="sample.id")

	# new outfile name
	phensavenm1 <- basename(phenfl)
	phensavenm2 <- gsub(".csv", "", phensavenm1)
	phensavenm3 <- paste0(phensavenm2, "_w_", gn, "_", to_process, ".csv")
	#phensavenm3 <- paste0(phensavenm2, "_w_", gn, "_noncoding.csv")
	phensavenm4 <- paste0(savedir, phensavenm3)

	message("Saving ", phensavenm4)
	write.csv(phen_w_vntcols, phensavenm4, row.names=F)
}



#vntcols_use <- vntcols[, grep(paste(genes_include, collapse="|"), names(vntcols))]
#vntcols_use$sample.id <- rownames(vntcols_use)
#
#phen_w_vntcols <- phen %>% left_join(vntcols_use, by="sample.id")
#
## new outfile name
#phensavenm1 <- basename(phenfl)
#phensavenm2 <- gsub(".csv", "", phensavenm1)
#phensavenm3 <- paste0(phensavenm2, "_w_addtl_noncod_vnt_cols.csv")
#phensavenm4 <- paste0("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/", phensavenm3)
#
#message("Saving ", phensavenm4)
#write.csv(phen_w_vntcols, phensavenm4, row.names=F)

