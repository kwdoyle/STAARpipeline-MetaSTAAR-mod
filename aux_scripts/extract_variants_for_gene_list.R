extractGenes <- function(x, genes) {
	dat <- get(load(x))
	dat_extract <- dat[genes]
	if (all(is.na(names(dat_extract)))) {
		# return something other than null if want to concat output list after?
		return(NULL)
	} else {
		message("Found genes")
		return(dat_extract) 
	}
}

#dir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Coding_Analysis/"
#dir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Coding_Analysis/"
# TODO extract these genes from the NONCODING output
#dir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/Discovery/Gene_Centric_Noncoding_Analysis/"
dir <- "/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/TOPMed/Gene_Centric_Noncoding_Analysis/"

genes_to_extract_df <- read.csv("/efs/garcia/users/kd2630/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/16_coding_hits_w_noncoding.csv")
genes <- unique(genes_to_extract_df$Gene.name)
# For coding genes, add C1QB to this
#genes <- c(genes, "C1QB")
# For noncoding, add C1QA and C1QC
genes <- c(genes, "C1QA", "C1QC")

sumstatfls <- list.files(dir, pattern="*sumstat*", full.names=T)
savenm <- gsub("_sumstat_[0-9]+\\.Rdata$", "", basename(sumstatfls[1]))
#savenm <- paste0(savenm, "_sumstat_for_extracted_genes_update.Rdata")
#savenm <- paste0(savenm, "_sumstat_for_16_coding_w_noncoding_genes.Rdata")
#savenm <- paste0(savenm, "_sumstat_for_16_coding_w_noncoding_genes_w_C1QB.Rdata")
savenm <- paste0(savenm, "_sumstat_for_16_coding_w_noncoding_genes_w_C1QA_C1QC.Rdata")
savedir <- paste0(dir, "/extracted_gene_output/")
dir.create(savedir, showWarnings=F)
savenm_full <- paste0(savedir,  savenm)
message("Will save output: ", savenm_full)

# per file, check if gene(s) to extract are in it and, if so, store output
#genes <- c("SAMD9L", "PPP1R3B")
#genes <- unique(c('PPP1R3B', 'TERT', 'SPDL1', 'SAMD9L', 'GOT2', 'MUC5B', 'CDCA2', 'ARL2', 'RCC1L', 'RCC1L', 'RCC1L', 'CDCA3',
#	   'ARL2-SNX15', 'TERT', 'TERT', 'P3H3', 'CNRIP1', 'C1QB',
#	  # add these additional genes
#	   'C1QC', 'C1QA'))
message("Genes to extract: ", paste(genes, collapse=", "))


sumstat_genes <- lapply(sumstatfls, extractGenes, genes=genes)
# get non-null elements
sumstat_genes_filt <- Filter(Negate(is.null), sumstat_genes)
sumstat_genes_filt <- unlist(sumstat_genes_filt, recursive=F)
# may be NA things in-between found genes--remove them again
sumstat_genes_filt2 <- Filter(Negate(is.null), sumstat_genes_filt)

message("Genes found: ", paste(names(sumstat_genes_filt2), collapse=", "))

save(sumstat_genes_filt2, file=savenm_full)

