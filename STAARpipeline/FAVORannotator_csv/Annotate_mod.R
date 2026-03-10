rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### DB split information 
#file_DBsplit <- "/home/STAAR/STAARpipeline-Tutorial/FAVORannotator_csv/FAVORdatabase_chrsplit.csv"

### anno channel (subset)
# this can stay hard-coded for now..
#anno_colnum <- c(1,8:12,15,16,19,23,25:36)
# Add in index of alphamissense column
anno_colnum <- c(1,8:12,15,16,19,23,25:36,38)
# 190 should be the last column in the data file
# !! Can now use hard-coded indices above, as they match with these new favor db downloads from the website

chr <- as.numeric(commandArgs(TRUE)[1])

xsv <- commandArgs(TRUE)[2]
output_path <- commandArgs(TRUE)[3]
DB_path <- commandArgs(TRUE)[4]
basedir <- commandArgs(TRUE)[5]

print(paste("xsv:", xsv))
print(paste("save path:", output_path))
print(paste("FAVOR db path:", DB_path))
print(paste("chromosome:", chr))

### DB split information 
file_DBsplit <- paste0(basedir, "/FAVORannotator_csv/FAVORdatabase_chrsplit.csv")
custom_annotation_path <- "~/STAAR_10k/10k_Cohort_codingonly/"

###########################################################################
#           Main Function 
###########################################################################

### chromosome number
## annotate (seperate)
DB_info <- read.csv(file_DBsplit,header=TRUE)
chr_splitnum <- sum(DB_info$Chr==chr)


for(kk in 1:chr_splitnum)
{
	print(kk)

	#system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Anno_chr",chr,"_",kk,".csv"))
	# first join anno from FAVOR
	system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/VarInfo_chr",chr,"_",kk,".csv variant_vcf ",DB_path,"/chr",chr,"_",kk,".csv > ",output_path,"chr",chr,"/Temp_FAVOR_chr",chr,"_",kk,".csv"))

	# then join custom (alphamissense) anno
	#system(paste0(xsv," join --left VarInfo chromosome position ref_vcf alt_vcf ",output_path,"chr",chr,"/Temp_FAVOR_chr",chr,"_",kk,".csv variant_vcf chromosome position ref_vcf alt_vcf ",custom_annotation_path,"/chr",chr,"_alphamissense.csv", " > ",output_path,"chr",chr,"/Anno_chr",chr,"_",kk,".csv"))  
	system(paste0(xsv," join --left VarInfo ",output_path,"chr",chr,"/Temp_FAVOR_chr",chr,"_",kk,".csv variant_vcf ",custom_annotation_path,"/chr",chr,"_alphamissense.csv", " > ",output_path,"chr",chr,"/Anno_chr",chr,"_",kk,".csv"))  
}

## merge info
Anno <- paste0(output_path,"chr",chr,"/Anno_chr",chr,"_",seq(1:chr_splitnum),".csv ")
merge_command <- paste0(xsv," cat rows ",Anno[1])

for(kk in 2:chr_splitnum)
{
	merge_command <- paste0(merge_command,Anno[kk])
}

merge_command <- paste0(merge_command,"> ",output_path,"chr",chr,"/Anno_chr",chr,".csv")

system(merge_command)

## subset
anno_colnum_xsv <- c()
for(kk in 1:(length(anno_colnum)-1))
{
	anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
}
anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])

system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,"chr",chr,"/Anno_chr",chr,".csv > ",output_path,"chr",chr,"/Anno_chr",chr,"_STAARpipeline.csv"))

