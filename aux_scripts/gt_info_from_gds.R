library(SeqArray)

#vntfl <- "~/MetaSTAAR_discovery_topmed_noncoding_gene_hit_variants_onlysig.xlsx"
vntfl <- "~/noncoding_variants_for_the_16_genes.xlsx"


#basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/10k_Cohort_rm_telo_qv/"
##savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Coding_variant_gt_matrices/"
#savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_variant_gt_matrices/"
#nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel_All_rm_qv.Rdata")  # To get sample IDs

basedir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/"
#savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Coding_variant_gt_matrices/"
savedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_variant_gt_matrices/"
nullmodelfl <- paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata") # To get sample IDs

message("Using:")
message("basedir: ", basedir)
message("savedir: ", savedir)
message("null model: ", nullmodelfl)

dir.create(savedir, showWarnings = F, recursive = T)


# TODO do this for the noncoding variants too, to condition the gene centric coding hits on first
#vnts <- readxl::read_xlsx("~/MetaSTAAR_discovery_topmed_coding_gene_hit_variants_onlysig.xlsx")
vnts <- readxl::read_xlsx(vntfl)

agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
obj_nullmodel <- get(load(paste0(nullmodelfl))) 

samp_ids <- obj_nullmodel$id_include

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
	# get GT matrix
	gts <- seqGetData(genofile, "genotype")
	gts_summed <- apply(gts, 2:3, sum, na.rm=T)
	rownames(gts_summed) <- samp_ids
	colnames(gts_summed) <- vnt_ids_use

	# get favor annotations
	#anno_nms <- ls.gdsn(index.gdsn(genofile, "annotation/info/FunctionalAnnotation/"))
	#annos.l <- lapply(anno_nms, function(x) seqGetData(genofile, paste0("annotation/info/FunctionalAnnotation/", x)))
	#names(annos.l) <- anno_nms

	save(gts_summed, file=paste0(savedir, "/chr", chr, "_variant_gt_matrix.rds"))

	seqClose(genofile)
}

