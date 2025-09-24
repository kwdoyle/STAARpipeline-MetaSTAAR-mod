rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
#file_DBsplit <- "/home/STAAR/STAARpipeline-Tutorial/FAVORannotator_csv/FAVORdatabase_chrsplit.csv"

chr <- as.numeric(commandArgs(TRUE)[1])

dir_geno <- commandArgs(TRUE)[2]
gds_file_name_1 <- commandArgs(TRUE)[3]
gds_file_name_2 <- commandArgs(TRUE)[4]
output_path <- commandArgs(TRUE)[5]
basedir <- commandArgs(TRUE)[6]

print(paste("base dir:", basedir))
print(paste("gds dir:", dir_geno))
print(paste("gds file name 1:", gds_file_name_1))
print(paste("gds file name 2:", gds_file_name_2))
print(paste("save path:", output_path))
print(paste("chromosome:", chr))

file_DBsplit <- paste0(basedir, "/FAVORannotator_csv/FAVORdatabase_chrsplit.csv")

###########################################################################
#           Main Function 
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
position <- as.integer(seqGetData(genofile, "position"))
REF <- as.character(seqGetData(genofile, "$ref"))
ALT <- as.character(seqGetData(genofile, "$alt"))

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
	print(kk)

	VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
	VarInfo <- data.frame(VarInfo)
	
	write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
}

