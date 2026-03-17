#!/bin/bash

#vcfnm=anno.gencode.coding.curegn.wgs.freeze.2a.chr22.hg38.filtered.samp.list.pruned.vcf.gz
for chr in {1..22}; do
	echo $chr

	vcfdir=/mnt/cuimc-med-pulm-general/vcfs/10k_hg38/samples_to_keep/GENCODE_coding_subset_100bp_exon_flank/annotated_vep115/
	vcfnm=anno.gencode.coding.curegn.wgs.freeze.2a.chr${chr}.hg38.filtered.samp.list.pruned.vcf.gz
	vcfl=${vcfdir}/${vcfnm}

	#outdir=~/STAAR_10k/10k_Cohort_codingonly/
	outdir=/mnt/tmpstore/10k_Cohort_codingonly/

	# %AlphaMissense_pred %AlphaMissense_rankscore

	# Print header
	#echo "variant_vcf,chromosome,position,ref_vcf,alt_vcf,alphamissense" > ${outdir}/chr22_alphamissense.csv
	echo "variant_vcf,alphamissense" > ${outdir}/chr${chr}_alphamissense.csv
	#echo "variant_vcf,alphamissense" > ${outdir}/chr22_alphamissense.csv

	bcftools +split-vep \
		${vcfl} \
		-d \
		-s worst \
		-f '%CHROM-%POS-%REF-%ALT,%AlphaMissense_score' \
		-A tab \
		| awk -F',' '{
		    n = split($2, a, "&");
		    max = a[1];
		    for (i = 2; i <= n; i++) {
			if (a[i] > max) max = a[i];
		    }
		    print $1","max;
		}' >> ${outdir}/chr${chr}_alphamissense.csv
done


# don't need the separate other fields
       	#-f '%CHROM-%POS-%REF-%ALT,%CHROM,%POS,%REF,%ALT,%AlphaMissense_score' \
#bcftools +split-vep -l ${vcfl}

# note:
# using without -d gives ~24k lines
# while using with -d gives ~ 283k lines
# but using with -d AND '-s worst' gives ~24k lines

       	#-d \


       	#-f '%CHROM-%POS-%REF-%ALT,%CHROM,%POS,%REF,%ALT,%AlphaMissense_score' \
       	#-f '%CHROM-%POS-%REF-%ALT,%CHROM,%POS,%REF,%ALT,%AlphaMissense_score,%Feature' \

	# testing:
	#-f '%CHROM-%POS-%REF-%ALT,%Feature,%ENSP,%HGVSp,%Protein_position,%Amino_acids,%Consequence,%BIOTYPE,%MANE_SELECT,%CANONICAL,%AlphaMissense_score' \
	#-f '%CHROM-%POS-%REF-%ALT,%SYMBOL,%Feature,%ENSP,%SWISSPROT,%TREMBL,%UNIPARC,%HGVSp,%Protein_position,%Amino_acids,%Consequence,%AlphaMissense_score' \
	#-f '%CHROM-%POS-%REF-%ALT,%SYMBOL,%Feature,%ENSP,%UNIPROT_ISOFORM,%SWISSPROT,%TREMBL,%UNIPARC,%HGVSp,%Protein_position,%Amino_acids,%Consequence,%AlphaMissense_pred,%AlphaMissense_rankscore,%AlphaMissense_score' \
	#-f '%CHROM-%POS-%REF-%ALT,%SYMBOL,%Feature,%ENSP,%UNIPROT_ISOFORM,%SWISSPROT,%TREMBL,%UNIPARC,%HGVSp,%Protein_position,%Amino_acids,%Consequence,%AlphaMissense_score' \
       	#-A tab | grep "22-15528267-A-G" 
#grep -e "22-15528267-A-G" -e "22-15528287-C-G" -e "22-15528306-T-A" -e "22-15528309-G-A" -e "22-15528315-A-T" -e "22-15528316-C-G" -e "22-15528331-G-T" -e "22-15528345-G-T" -e "22-15528351-T-C" -e "22-15528543-C-T"
