
dircache=/mnt/cuimc-med-pulm-general/vepcache114/

~/ensembl-vep/vep -i $file \
	-o ${savedir}/${saveflnm} \
	--tab \
	--cache --offline \
	--dir_cache ${dircache} \
	--assembly GRCh38 \
	-e --no_stats --buffer_size 5000

#	--vcf --force_overwrite \
	#--compress_output bgzip \
