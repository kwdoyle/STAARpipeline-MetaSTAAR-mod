library(SeqArray)

#basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/10k_Cohort_rm_telo_qv/"
#savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Coding_variant_gt_matrices/"
#nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel_All_rm_qv.Rdata")

basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/"
savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Coding_variant_gt_matrices/"
nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata")

message("Using:")
message("basedir: ", basedir)
message("savedir: ", savedir)
message("null model: ", nullmodelfl)


dir.create(savedir, showWarnings = F, recursive = T)
coding_vnts <- readxl::read_xlsx("~/MetaSTAAR_discovery_topmed_coding_gene_hit_variants_onlysig.xlsx")

agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

samp_ids <- obj_nullmodel$id_include

# make pasted ref,alt column to match allele format in GDS
coding_vnts$alleles <- paste(coding_vnts$ref, coding_vnts$alt, sep=",")

# test extracting gts from one chr gds file
chrs_use <- unique(coding_vnts$chr)

cat("\n")
for (chr in chrs_use) {
	message("Chromosome ", chr)
	# NOTE: I assume that, when using the full list of variants that is between disc and tm,
	# the subsetted list of variants I get from the GDS using 'vnts_to_extract' should
	# match the number I see when I filter out those where MAC_cohortX (for disc or tm) is 0
	vnts_to_extract_cols <- coding_vnts[coding_vnts$chr == chr, c("chr","pos","alleles")]
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
	# get GT matrix
	gts <- seqGetData(genofile, "genotype")
	gts_summed <- apply(gts, 2:3, sum, na.rm=T)
	rownames(gts_summed) <- samp_ids
	colnames(gts_summed) <- vnt_ids_use

	save(gts_summed, file=paste0(savedir, "/chr", chr, "_variant_gt_matrix.rds"))

	seqClose(genofile)
}

