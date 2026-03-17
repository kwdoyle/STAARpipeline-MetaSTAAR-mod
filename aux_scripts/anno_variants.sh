#!/bin/bash

chr=$1
type=topmed  #topmed  #discovery

# activate vep conda environment first.
# this was previously set up to run the loftee plugin:

# conda activate vep-loftee

#vepuse=~/ensembl-vep/vep
vepuse=/home/kd2630/miniconda3/envs/vep-loftee/bin/vep

# still use the vepcache stored on the S3, since this is too large to copy to tmpstore
dircache=/mnt/cuimc-med-pulm-general/vepcache115/

if [[ $type == "discovery" ]]; then
	vcfdir=/mnt/cuimc-med-pulm-general/vcfs/10k_hg38/samples_to_keep/GENCODE_coding_subset_100bp_exon_flank/
	#vcfdir=/home/kd2630/tmp/
	vcf1=gencode.coding.curegn.wgs.freeze.2a.chr
	vcf2=.hg38.filtered.samp.list.pruned.vcf.gz
elif [[ $type == "topmed" ]]; then
	vcfdir=/mnt/cuimc-med-pulm-general/vcfs/topmed_full_concat/GENCODE_coding_subset_100bp_exon_flank/
	vcf1=gencode.coding.topmed_chr
	vcf2=.vcf.gz
fi

mainsavedir=${vcfdir}annotated_vep115/
savedir=/mnt/tmpstore/annotated_vep115/
mkdir $mainsavedir
mkdir $savedir



file=${vcf1}${chr}${vcf2}
echo $file
saveflnm=anno."${file##*/}"

# temporarily copy the vcf to /mnt/tmpstore/ and operate on it there
s3_vcfdir=${vcfdir/\/mnt\//s3://}
tmpfile=/mnt/tmpstore/${file}
aws s3 cp ${s3_vcfdir}${file} ${tmpfile}


# dbNSFP database is restoring from deep archive. will be available until 3/14/26.
# has been copied over to /mnt/tmpstore anyway, so its available to use until I delete it.
#$vepuse -i ${vcfdir}/${file} \
$vepuse -i ${tmpfile} \
	-o ${savedir}/${saveflnm} \
	--vcf --force_overwrite \
	--cache --offline \
	--dir_cache ${dircache} \
	--assembly GRCh38 \
	--compress_output bgzip \
	-e --no_stats --buffer_size 5000 \
	--plugin dbNSFP,/mnt/tmpstore/dbNSFP4.7a_grch38.gz,CADD_raw,CADD_raw_rankscore,CADD_phred,REVEL_score,REVEL_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred

# mv to s3 when done
s3_mainsavedir=${mainsavedir/\/mnt\//s3://}
aws s3 mv ${savedir}${saveflnm} ${s3_mainsavedir}${saveflnm}











	#--plugin dbNSFP,/mnt/cuimc-med-pulm-general/dbNSFP/dbNSFP4.7a_grch38.gz,CADD_raw,CADD_raw_rankscore,CADD_phred,REVEL_score,REVEL_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred

	#--plugin LoF,loftee_path:/home/kd2630/loftee-1.0.4_GRCh38/,human_ancestor_fa:${loftee_ancest_fl},gerp_bigwig:${gerpfl} \
       	#--dir_plugins /home/kd2630/loftee-1.0.4_GRCh38/ \

#--buffer_size 5000
#~/ensembl-vep/vep -i $file \
#done

#for file in "$vcfdir"/*.vcf.gz; do

# filter vcf for variants with MAF < 0.01 before annotating
#	bcftools +fill-tags $file -- -t MAF |
#		bcftools view -i 'INFO/MAF<0.01' -Ov |
#~/ensembl-vep/vep -i /dev/stdin \

#file=/home/kd2630/STAAR_10k/MetaSTAAR/coding_gene_hit_variants_onlysig3.tsv
#file=/home/kd2630/STAAR_10k/MetaSTAAR/noncoding_gene_hit_variants_onlysig3.tsv
#saveflnm=coding_gene_hit_variant_anno.tsv
#saveflnm=noncoding_gene_hit_variant_anno.tsv
#fastafl=/mnt/cuimc-med-pulm-general/vepcache114/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#loftee_ancest_fl=/mnt/cuimc-med-pulm-general/LOFTEE/human_ancestor.fa.gz
#gerpfl=/mnt/cuimc-med-pulm-general/LOFTEE/gerp_conservation_scores.homo_sapiens.GRCh38.bw

