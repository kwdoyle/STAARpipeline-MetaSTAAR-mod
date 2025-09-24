###########################################################
# fit STAAR null model
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 01/06/2023
###########################################################
rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Phenotype file
#phenotype <- read.csv("/home/STAAR/TOPMed_Full_Cohort/topmed_input_model_data_cleaned.csv")
## (sparse) GRM file
#sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))
## file directory for the output file 
#output_path <- "/home/STAAR/TOPMed_Full_Cohort/staar_null_model/"
## output file name
#output_name <- "obj_nullmodel.Rdata"

phenotype_fl <- commandArgs(TRUE)[1]
output_name <- commandArgs(TRUE)[2]
sgrm_fl <- commandArgs(TRUE)[3]

phenotype <- read.csv(phenotype_fl)
sgrm <- get(load(sgrm_fl))

print(paste("Using file", phenotype_fl))
print(paste("Using GRM file", sgrm_fl))
print(paste("Saving to", output_name))

#dir.create(output_path)

# Need to remove the "0_" from sample IDs in GRM
rownames(sgrm) <- gsub("0_", "", rownames(sgrm))
colnames(sgrm) <- gsub("0_", "", colnames(sgrm))

###########################################################
#           Main Function 
###########################################################
## fit null model
#obj_nullmodel <- fit_nullmodel(is_IPF~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
#                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
#                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)

obj_nullmodel <- fit_nullmodel(is_IPF~Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,  #+as.factor(rs35705950), 
                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)



save(obj_nullmodel,file=output_name)

