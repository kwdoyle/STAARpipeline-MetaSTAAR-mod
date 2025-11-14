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
# If two covariates are virtually the same, remove one of them
detect_colinear_vars <- function(data, threshold = 0.99) {  
  # Select only numeric/binary columns  
  numeric_cols <- sapply(data, function(x) is.numeric(x) || all(x %in% c(0, 1, NA)))  
  data_numeric <- data[, numeric_cols, drop = FALSE]  
    
  # Calculate correlation matrix  
  cor_matrix <- cor(data_numeric, use = "pairwise.complete.obs")  
    
  # Find pairs with correlation above threshold  
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- 0  
  high_cor <- which(abs(cor_matrix) > threshold, arr.ind = TRUE)  
    
  if (nrow(high_cor) == 0) {  
    return(list(to_remove = character(0), pairs = data.frame()))  
  }  
    
  # Create pairs dataframe  
  pairs_df <- data.frame(  
    var1 = rownames(cor_matrix)[high_cor[, 1]],  
    var2 = colnames(cor_matrix)[high_cor[, 2]],  
    correlation = cor_matrix[high_cor],  
    stringsAsFactors = FALSE  
  )  
    
  # Identify variables to remove (keep first occurrence)  
  to_remove <- unique(pairs_df$var2)  
    
  return(list(to_remove = to_remove, pairs = pairs_df))  
}  

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

# New: check if any variables are all the same for all samples. If so, remove that variable.
same_check <- unlist(lapply(phenotype, function(x) length(unique(x)) == 1))
if (any(same_check)) {
	rm_col <- names(which(same_check))
	message("All values identical for ", rm_col, " -- omitting from model")
	phenotype <- phenotype[, -which(names(phenotype) %in% rm_col)]
}

# check for colinear terms
colin_chk <- detect_colinear_vars(phenotype, threshold=0.95)
if (length(colin_chk$to_remove) > 0) {  
	message("Colinear variables found; removing ", colin_chk$to_remove)
	phenotype <- phenotype[, -which(names(phenotype) %in% colin_chk$to_remove)]
}

# new: dynamically create model formula using whichever 'has_vnt_in' columns are present
fixed_covariates <- c("Age", "Sex", paste0("PC", 1:10))
vnt_cols <- grep("^has_vnt_in_", names(phenotype), value = TRUE)
# combine
predictors <- c(fixed_covariates, vnt_cols)

form <- as.formula(paste("is_IPF ~", paste(predictors, collapse = " + ")))

obj_nullmodel <- fit_nullmodel(form, 
                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)


# old:
#obj_nullmodel <- fit_nullmodel(is_IPF~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, 
#                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
#                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)

# test
#obj_nullmodel <- fit_nullmodel(is_IPF~Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + has_vnt_in_GOT2_promoter_CAGE + has_vnt_in_GOT2_enhancer_CAGE, 
#                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
#                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)

# previous latest way:
#obj_nullmodel <- fit_nullmodel(is_IPF~Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,  #+as.factor(rs35705950), 
#                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
#                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)

# adjusting gene-centric coding output on presense of any noncoding variant in same gene
#obj_nullmodel <- fit_nullmodel(is_IPF~Age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ 
#			       has_vnt_in_ARL2_promoter_CAGE + has_vnt_in_GOT2_promoter_CAGE + has_vnt_in_GOT2_enhancer_CAGE +
#			       has_vnt_in_P3H3_promoter_DHS + has_vnt_in_SAMD9L_upstream + has_vnt_in_SAMD9L_promoter_DHS + has_vnt_in_SAMD9L_enhancer_DHS,
#                               data=phenotype, kins=sgrm, use_sparse=TRUE, kins_cutoff=0.022, id="sample.id",
#                               groups=NULL, family=binomial(link="logit"), verbose=TRUE)



save(obj_nullmodel,file=output_name)

