library(SeqArray)

cohort <- "discovery"  # "discovery"  "topmed"

#vntfl <- "~/coding_variants_for_the_18_genes_w_c1q.xlsx"
#vntfl <- "~/coding_variants_for_the_16_genes.xlsx"
#vntfl <- "~/noncoding_variants_for_the_16_genes.xlsx"
vntfl <- "~/noncoding_variants_for_the_18_genes_w_c1q.xlsx"
#process <- "coding"  # "coding"  "noncoding"

savedirnm <- "Gene_Centric_Noncoding_sighits_variant_gt_matrices"  #"Gene_Centric_Coding_sighits_variant_gt_matrices" # "Gene_Centric_Noncoding_sighits_variant_gt_matrices" #"Gene_Centric_Noncoding_variant_gt_matrices"

# don't use these files:
#vntfl <- "~/MetaSTAAR_discovery_topmed_coding_gene_hit_variants_onlysig.xlsx"
#vntfl <- "~/MetaSTAAR_discovery_topmed_noncoding_gene_hit_variants_onlysig.xlsx"

# ---

# TODO extract variants for all noncoding and coding signig staar hits
# using the tables in the 'send_to_gareth' folder
# also set this up to set the below params better for the cohorts and for coding or noncoding
basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/"

# set directories and files
if (cohort == "discovery") {
	subdir <- "/10k_Cohort_rm_telo_qv/"
	grmpth <- paste0(subdir, "/staar_null_model/obj_nullmodel_All_rm_qv.Rdata")
	phenpth <- paste0(subdir, "/discovery_cohort_844_controls_rm_telo_gene_carriers_input_model_data.csv")
	cohortdir <- "Discovery"
} else if (cohort == "topmed") {
	subdir <- "/TOPMed_Full_Cohort_grm/"
	grmpth <- paste0(subdir, "/staar_null_model/obj_nullmodel.Rdata")
	phenpth <- paste0(subdir, "/topmed_input_model_data_rv_rm.csv")
	cohortdir <- "TOPMed"
} else {
	stop("cohort incorrectly specified.")
}

nullmodelfl <- paste0(basedir, grmpth) 
savedir <- paste0(basedir, "/MetaSTAAR_Discovery_TOPMed/", cohortdir, "/",  savedirnm, "_", cohortdir)
phenfl <- paste0(basedir, phenpth)

#basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/10k_Cohort_rm_telo_qv/"
#savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Coding_variant_gt_matrices_Discovery_Update/"
##savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_variant_gt_matrices_Discovery_Update/"
#nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel_All_rm_qv.Rdata")  # To get sample IDs
#phenfl <- paste0(basedir, "/discovery_cohort_844_controls_rm_telo_gene_carriers_input_model_data.csv")

#basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/"
##savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Coding_variant_gt_matrices_TOPMed_Update/"
#savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_variant_gt_matrices_TOPMed_Update/"
#nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata") # To get sample IDs
#phenfl <- paste0(basedir, "/topmed_input_model_data_rv_rm.csv")

message("Using:")
#message("basedir: ", basedir)
message("savedir: ", savedir)
message("null model: ", nullmodelfl)
message("variant file: ", vntfl)
AFdir <- paste0(savedir, "/case_ctrl_AF_dfs/")
message("Will save case/control AF data frames to: ", AFdir)

dir.create(savedir, showWarnings = F, recursive = T)
dir.create(AFdir, showWarnings = F, recursive = T)


# TODO do this for the noncoding variants too, to condition the gene centric coding hits on first
#vnts <- readxl::read_xlsx("~/MetaSTAAR_discovery_topmed_coding_gene_hit_variants_onlysig.xlsx")
vnts <- readxl::read_xlsx(vntfl)

agds_dir <- get(load(paste0(basedir, subdir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

samp_ids <- obj_nullmodel$id_include
# TODO need to use actual case/control list for TOPMed
pheno <- read.csv(phenfl)
cases <- pheno$sample.id[pheno$is_IPF == 1]
controls <- pheno$sample.id[pheno$is_IPF == 0]
#cases <- grep("IPF|Pul", samp_ids, value=T)
#controls <- grep("IPF|Pul", samp_ids, value=T, invert=T) 

# make pasted ref,alt column to match allele format in GDS
vnts$alleles <- paste(vnts$ref, vnts$alt, sep=",")

# test extracting gts from one chr gds file
chrs_use <- unique(vnts$chr)

# TODO also extract the favor annotations

cat("\n")
for (chr in chrs_use) {
	message("Chromosome ", chr)
	# NOTE: I assume that, when using the full list of variants that is between disc and tm,
	# the subsetted list of variants I get from the GDS using 'vnts_to_extract' should
	# match the number I see when I filter out those where MAC_cohortX (for disc or tm) is 0
	vnts_to_extract_cols <- vnts[vnts$chr == chr, c("chr","pos","alleles")]
	vnts_to_extract <- paste(vnts_to_extract_cols$chr, vnts_to_extract_cols$pos, vnts_to_extract_cols$alleles, sep="-")

	#chr = chrs_use[1]
	agds.path <- agds_dir[chr]
	genofile <- seqOpen(agds.path)

	# check variant.ids
	vids <- seqGetData(genofile, "variant.id")  # these are just numbers--useless
	chrs <- seqGetData(genofile, "chromosome")
	poss <- seqGetData(genofile, "position")
	alleles <- seqGetData(genofile, "allele")
	# make actual variant id
	gds_vnt_ids <- paste(chrs, poss, alleles, sep="-")
	# check sample ids
	#samps <- seqGetData(genofile, "sample.id")

	vnt_idx_use <- which(gds_vnt_ids %in% vnts_to_extract)
	# use these to label the gt matrix
	vnt_ids_use <- gds_vnt_ids[vnt_idx_use]
	# use the corresponding variant id 'number' that corresponds to the chr pos ref and alt to filter the variants
	vids_use <- vids[vnt_idx_use]


	# Filter gds for samples and variants
	seqSetFilter(genofile, variant.id = vids_use, sample.id = samp_ids)
	# remake sample IDs to match order in gds
	actual_samps <- seqGetData(genofile, "sample.id") 
	# remake variant IDs
	# NOTE: variant IDs match, since they were found by finding the index in the gds in the first place,
	# while sample IDs were not
#	vids2 <- seqGetData(genofile, "variant.id")  # these are just numbers--useless
#	chrs2 <- seqGetData(genofile, "chromosome")
#	poss2 <- seqGetData(genofile, "position")
#	alleles2 <- seqGetData(genofile, "allele")
#	# make actual variant id
#	gds_vnt_ids2 <- paste(chrs2, poss2, alleles2, sep="-")
	# check sample ids
	#samps <- seqGetData(genofile, "sample.id")

#	vnt_idx_use2 <- which(gds_vnt_ids2 %in% vnts_to_extract2)
#	# use these to label the gt matrix
#	vnt_ids_use2 <- gds_vnt_ids2[vnt_idx_use2]

	# get GT matrix
	gts <- seqGetData(genofile, "genotype")
	gts_summed <- apply(gts, 2:3, sum, na.rm=T)
	rownames(gts_summed) <- actual_samps  #samp_ids
	colnames(gts_summed) <- vnt_ids_use

	# get favor annotations
	#anno_nms <- ls.gdsn(index.gdsn(genofile, "annotation/info/FunctionalAnnotation/"))
	#annos.l <- lapply(anno_nms, function(x) seqGetData(genofile, paste0("annotation/info/FunctionalAnnotation/", x)))
	#names(annos.l) <- anno_nms

	# temporarily do not save when re-saving the AF data frames below
	save(gts_summed, file=paste0(savedir, "/chr", chr, "_variant_gt_matrix.rds"))

	# get case and control ACs, AFs
	samps <- seqGetData(genofile, "sample.id")
	case_idx <- which(samps %in% cases)  #match(cases, samps)
	control_idx <- which(samps %in% controls) #match(controls, samps)

	# summed gts accounting for missingness
	gts_summed2 <- apply(gts, c(2,3), function(a) {
		  if (any(is.na(a))) NA_integer_ else sum(a)
		})
	rownames(gts_summed2) <- samps
	colnames(gts_summed2) <- vnt_ids_use

	called_alleles <- apply(gts, c(2,3), function(a) sum(!is.na(a)))

	gts_summed2_case <- gts_summed2[case_idx, , drop=F]
	gts_summed2_control <- gts_summed2[control_idx, , drop=F]
	called_allele_case <- called_alleles[case_idx, , drop=F]
	called_allele_control <- called_alleles[control_idx, , drop=F]

	AC_case <- colSums(gts_summed2_case, na.rm=T)
	AN_case <- colSums(called_allele_case)
	AF_case <- AC_case / pmax(AN_case, 1)

	AC_control <- colSums(gts_summed2_control, na.rm=T)
	AN_control <- colSums(called_allele_control)
	AF_control <- AC_control / pmax(AN_control, 1)

	# make AF data frame
	vntAF_df <- data.frame(variantID = vnt_ids_use,
				AC_case=AC_case, AN_case=AN_case, AF_case=AF_case,
				AC_control=AC_control, AN_control=AN_control, AF_control=AF_control)

	write.csv(vntAF_df, file=paste0(AFdir, "/chr", chr, "_case_ctrl_AFs.csv"), row.names=F)

	seqClose(genofile)
}

