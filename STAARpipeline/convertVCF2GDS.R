# Takes a VCF file as input and outputs it to GDS format.
# A wrapper for the seqVCF2GDS function in the SeqArray package.
# options(repos = c(CRAN = "http://cran.revolutionanalytics.com"))
# VCF Input: can be one VCF or several, as a character vector
# Author(s): Michael R. Brown, Jennifer A. Brody


args <- commandArgs(trailingOnly = TRUE)
print(args)
import.format <- args[1] # characters, the variable name(s) in the FORMAT field for import; or NULL for all variables
fileformat <- args[2] # either 'vcf' or 'bcf' (or 'bed')
base.filename <- args[3] # output GDS file name
num.files <- as.numeric(args[4]) # number of input VCF files
vcf.file <- args[5:(5+num.files-1)] # input VCF file names
savedir <- args[6]
print(vcf.file)
print(paste("Saving to", savedir))

# GDS Ouput
gds.file <- paste0(base.filename,".gds")
# Make parent directory if it doesn't exist yet
dir.create(dirname(gds.file), showWarnings=F)
gds.file.basename <- basename(gds.file)

# check if gds file already exists and quit if it does (can later allow to overwrite file if specify)
if (file.exists(gds.file)) {
	stop("GDS file already created. Exiting..")
}

cmd <- paste0("aws s3 mv ", gds.file, " s3://cuimc-med-pulm-general/STAAR/", savedir, "/",  gds.file.basename)
cmd_clean <- gsub("(?<!:)//+", "/", cmd, perl = TRUE)

print(paste("Will run command:", cmd_clean))

library(SeqArray)

# Conversion -- This could take some time (if the input files are large)
# note: on the sge cluster, using the parallel argument causes the job to instantly crash.
#nSlots <- system("nproc",intern=T)
#nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
#if (nThreads == 0) nThreads <- 1
#message(paste("Running with", nThreads,"thread(s)."))
#Sys.setenv(MKL_NUM_THREADS=nThreads)


if (fileformat == 'bcf') {
        ## is this a bcf file?
        ## use bcftools to read text
        message("converting BCF")
        message("BCF always single thread")
	seqBCF2GDS(vcf.file,gds.file,fmt.import=import.format, storage.option="LZMA_RA",bcftools="bcftools")
} else if (fileformat == 'vcf') {
    message("converting VCF")
	seqVCF2GDS(vcf.file,gds.file,fmt.import=import.format, storage.option="LZMA_RA")  #,parallel=nThreads)
} else if (fileformat == 'bed') {
	message("converting BED")
	seqBED2GDS(vcf.file,gds.file,fmt.import=import.format, storage.option="LZMA_RA")  #,parallel=nThreads)
} else {
	stop("file format must be 'bcf', 'vcf', or 'bed'")
}
# Quick summary for log's sake
g <- seqOpen(gds.file)
seqSummary(g)

system(cmd_clean) 
#system(paste0("aws s3 mv ", gds.file, " s3://cuimc-med-pulm-general/STAAR/TOPMed_codingonly/", gds.file.basename)) 

