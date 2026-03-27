rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

# these will always be the same
anno_file_name_1 <- "Anno_chr"
anno_file_name_2 <- "_STAARpipeline.csv"

chr <- as.numeric(commandArgs(TRUE)[1])

dir_geno <- commandArgs(TRUE)[2]
gds_file_name_1 <- commandArgs(TRUE)[3]
gds_file_name_2 <- commandArgs(TRUE)[4]
dir_anno <- commandArgs(TRUE)[5]

message("aGDS for chr ", chr)

###########################################################################
#           Main Function
###########################################################################

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)

### read annotation data
# NOTE: read_csv reads in too many columns incorrectly as logical. Just use data.table::fread instead
message("Reading ", paste0(dir_anno,"chr",chr,"/",anno_file_name_1,chr,anno_file_name_2))
FunctionalAnnotation <- data.table::fread(paste0(dir_anno,"chr",chr,"/",anno_file_name_1,chr,anno_file_name_2))

# set alphamissense to numeric, if it exists
if ('alphamissense' %in% names(FunctionalAnnotation)) {
	FunctionalAnnotation$alphamissense <- as.numeric(FunctionalAnnotation$alphamissense)
}

dim(FunctionalAnnotation)
print(str(FunctionalAnnotation))


## rename colnames
# Update: with the newer FAVOR db downloads, these renames are now accurate

# if these names don't exist in FunctionalAnnotation, then they won't be renamed.
# this way, if by some chance the 2nd, 7th, and 9th columns aren't these variables,
# the wrong thing won't get renamed
rn_idx1 <- grep("apc_conservation_v2", names(FunctionalAnnotation))
rn_idx2 <- grep("apc_local_nucleotide_diversity_v3", names(FunctionalAnnotation))
rn_idx3 <- grep("apc_protein_function_v3", names(FunctionalAnnotation))

print(paste("Renaming", colnames(FunctionalAnnotation)[rn_idx1], "to apc_conversion"))
colnames(FunctionalAnnotation)[rn_idx1] <- "apc_conservation"

print(paste("Renaming", colnames(FunctionalAnnotation)[rn_idx2], "apc_local_nucleotide_diversity"))
colnames(FunctionalAnnotation)[rn_idx2] <- "apc_local_nucleotide_diversity"

print(paste("Renaming", colnames(FunctionalAnnotation)[rn_idx3], "apc_protein_function"))
colnames(FunctionalAnnotation)[rn_idx3] <- "apc_protein_function"


## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)

### !! Cannot annotate gds while its in the S3 bucket,
# so copy it to a tmp folder first, annotate it, then copy it back over

# copy to tmp
system(paste0("cp ", gds.path, " /mnt/tmpstore/"))
gds.path.local <- paste0("/mnt/tmpstore/", basename(gds.path))

#genofile <- seqOpen(gds.path, readonly = FALSE)
genofile <- seqOpen(gds.path.local, readonly = FALSE)

# check that the variant order in FunctionalAnnotation matches that in the gds file.
# if it doesn't, throw an error.
pos = seqGetData(genofile, "position")
refalt = seqGetData(genofile, "allele")
refalt2 = gsub(",", "-", refalt)
varid_gds = paste(chr, pos, refalt2, sep="-")
# subset func anno data for the exact variants in the gds because, for some reason, there are still missmatches.
#gds_vnt_match_idx <- which(FunctionalAnnotation$VarInfo %in% varid_gds)
#FunctionalAnnotation_match <- FunctionalAnnotation[gds_vnt_match_idx, ]
# need all matches of varid_gds in FuncAnno--even if it duplicates the FuncAnno row..
gds_vnt_match_idx <- match(varid_gds, FunctionalAnnotation$VarInfo)
FunctionalAnnotation_match <- FunctionalAnnotation[gds_vnt_match_idx, ]

stopifnot(all(FunctionalAnnotation_match$VarInfo == varid_gds))

Anno.folder <- index.gdsn(genofile, "annotation/info")
# add replace=T to overwrite old annotations?
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation_match, compress="LZMA_ra", closezip=TRUE, replace=TRUE)

seqClose(genofile)

# remove extra slashes in command so that aws copies to correct location
# (and doesn't make stupid '/' directoires WITHIN a directory)
cmd <- paste0("aws s3 mv ", gds.path.local, " ",  gsub("/mnt/", "s3://", gds.path))
cmd_clean <- gsub("(?<!:)//+", "/", cmd, perl = TRUE)
system(cmd_clean)

