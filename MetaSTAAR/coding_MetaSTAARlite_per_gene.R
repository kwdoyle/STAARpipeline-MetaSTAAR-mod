library(STAAR)
library(STAARpipeline)
library(MetaSTAAR)
library(MetaSTAARlite)

## Directories of the study-specific summary statistics file folders
file.dir <- c("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Coding_Analysis_noncod_adj_null_mod_Update_w_c1q/",
              "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Coding_Analysis_noncod_adj_null_mod_Update_w_c1q/")
file.prefix <- c("10k_Cohort_rm_telo_qv_Gene_Coding_Analysis","TOPMed_Full_Cohort_grm_Gene_Coding_Analysis")
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
output_path <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Coding_Meta_Analysis_noncod_adj_null_mod_Update_w_c1q/"
dir.create(output_path, recursive=T, showWarnings=F)
## output file name
output_file_name <- "codingMeta"


results_coding <- c()

sumstat.file.path <- paste0(file.dir, file.prefix, "_sumstat.Rdata")
cov.file.path <- paste0(file.dir, file.prefix, "_cov.Rdata")
coding_sumstat_list <- sapply(sumstat.file.path, function(x) mget(load(x)), simplify = TRUE)
coding_cov_list <- sapply(cov.file.path, function(x) mget(load(x)), simplify = TRUE)

# genes to analyze
genes_to_analyze <- unique(unlist(lapply(coding_sumstat_list, names)))

for (gn_nm in genes_to_analyze) {
	gn <- strsplit(gn_nm, "_on_")[[1]][1]
	second_part <- sub(".*_on_", "", gn_nm)
	if (gn != second_part) {
		savenm <- paste(gn, "on", second_part, sep="_")
	} else {
		savenm <- gn
	}

	chr <- genes_info[genes_info$hgnc_symbol == gn, "chromosome_name"]
	message("Gene: ", gn)
	message("Chr: ", chr)

	  coding_sumstat_gene_list <- lapply(sumstat.file.path, function(x) {
	    coding_sumstat_list[[paste0(x,".coding_sumstat")]][[savenm]]  # [[gn]]
	  })
	  coding_cov_gene_list <- lapply(cov.file.path, function(x) {
	    coding_cov_list[[paste0(x,".coding_cov")]][[savenm]]  # [[gn]]
	  })

	  results <- coding_MetaSTAARlite(chr=chr, gene_name=gn, genes=genes_info,
					  sample.sizes=sample.sizes, coding_sumstat_gene_list=coding_sumstat_gene_list,
					  coding_cov_gene_list=coding_cov_gene_list,
					  cov_maf_cutoff=cov_maf_cutoff,
					  rare_maf_cutoff=0.01, rv_num_cutoff=2,
					  check_qc_label=TRUE, variant_type=variant_type,
					  Use_annotation_weights=Use_annotation_weights, Annotation_name=Annotation_name)

	# put into list with name of current gene
	results_l <- list(results)
	names(results_l) <- savenm #gn

	results_coding <- append(results_coding, results_l)

}

save(results_coding, file=paste0(output_path, output_file_name, ".Rdata"))

