
# NOTE: this script is used to combine all the separate gene outputs from the 'worker' conditional script.
# Ideally, the separate genes could all just be saved as one big list from that script itself.

path <- commandArgs(trailingOnly=T)[1]
message("Aggregating files in ", path)


covcond_fls <- list.files(path, pattern="cov_cond.Rdata", full.names=T)
cov_fls <- list.files(path, pattern="cov.Rdata", full.names=T)
sumstat_fls <- list.files(path, pattern="sumstat.Rdata", full.names=T)

if (grepl("NotCond", covcond_fls[1])) {
	replace_nm <- "NotCond"
} else if (grepl("Cond", covcond_fls[1])) {
	replace_nm <- "Cond"
} else {
	stop("Cannot find conditional or not-conditional output files")
}

covcond_dat_l <- lapply(covcond_fls, function(x) get(load(x)))
cov_dat_l <- lapply(cov_fls, function(x) get(load(x)))
sumstat_dat_l <- lapply(sumstat_fls, function(x) get(load(x)))

# reduce to single list with named elements for each gene
covcond_dat <- unlist(covcond_dat_l, recursive=F)
cov_dat <- unlist(cov_dat_l, recursive=F)
sumstat_dat <- unlist(sumstat_dat_l, recursive=F)

# extract base file names and save
#covcond_nm <- paste0(sub("_Cond_.*", "_Cond_", basename(covcond_fls[1])), "cov_cond.Rdata")
#cov_nm <- paste0(sub("_Cond_.*", "_Cond_", basename(cov_fls[1])), "cov.Rdata")
#sumstat_nm <- paste0(sub("_Cond_.*", "_Cond_", basename(sumstat_fls[1])), "sumstat.Rdata")
covcond_nm <- paste0(sub(paste0("_", replace_nm, "_.*"), paste0("_", replace_nm, "_"), basename(covcond_fls[1])), "cov_cond.Rdata")
cov_nm <- paste0(sub(paste0("_", replace_nm, "_.*"), paste0("_", replace_nm, "_"), basename(cov_fls[1])), "cov.Rdata")
sumstat_nm <- paste0(sub(paste0("_", replace_nm, "_.*"), paste0("_", replace_nm, "_"), basename(sumstat_fls[1])), "sumstat.Rdata")

message("Saving ", paste0(path, "/", covcond_nm))
message("Saving ", paste0(path, "/", cov_nm))
message("Saving ", paste0(path, "/", sumstat_nm))

save(covcond_dat, file=paste0(path, "/", covcond_nm), compress = "xz")
save(cov_dat, file=paste0(path, "/", cov_nm), compress = "xz")
save(sumstat_dat, file=paste0(path, "/", sumstat_nm), compress = "xz")

