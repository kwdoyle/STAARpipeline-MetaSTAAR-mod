
basedir <- "~/noncoding_telo/STAAR/MetaSTAAR_Discovery_TOPMed/"
#diruse <-  "/noncoding_gene_vnt_null_models/"  # "/coding_gene_vnt_null_models/"  # "/noncoding_gene_vnt_null_models/"
diruse <-  "/noncoding_gene_vnt_null_models_w_c1q/"  # "/coding_gene_vnt_null_models/"  # "/noncoding_gene_vnt_null_models/"
thedir <- paste0(basedir, diruse)

savenm <- "noncoding_vnt_categs_used_to_adj_w_c1q" #"coding_vnt_categs_used_to_adj_w_c1q"  # "coding_vnt_categs_used_to_adj"  # "noncoding_vnt_categs_used_to_adj"

mods <- list.files(thedir, pattern = "^discovery.*Rdata$", full.names=T)

getParams <- function(fl) {
	nm <- get(load(fl))
	allparams <- names(nm$coefficients)
	params <- grep("has_vnt", allparams, value=T)

	return(params)
}

vnts_adj <- unique(lapply(mods, getParams))
# assign gene names to each list element
gene_names <- sapply(vnts_adj, function(x) {
  unique(sub(".*_in_([A-Z0-9]+).*", "\\1", x))
})
names(vnts_adj) <- gene_names

# turn to data frame
vnts_adj_df <- stack(vnts_adj)
vnts_adj_df <- vnts_adj_df[,c(2,1)]
names(vnts_adj_df) <- c("Gene", "adjusted_for")

#print(vnts_adj)
write.csv(vnts_adj_df, paste0(basedir, "/", savenm, ".csv"), row.names=F)

