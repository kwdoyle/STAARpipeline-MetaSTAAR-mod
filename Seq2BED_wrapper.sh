#!/bin/bash
source /efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/makeenv2.sh

# This script sources my env and actually runs the conversion script

chr=$1
basedir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/
#basedir=~/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/GRM/

# --- change these:
#gdsdir=~/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/
gdsdir=/efs/garcia/users/kd2630/noncoding_telo/STAAR/10k_Cohort_rm_telo_qv/
savedir=${gdsdir}/beds_full/

gdsfl=${gdsdir}/10k_chr_${chr}.gds
savenm=${savedir}/10k_chr_${chr}.out
#gdsfl=${gdsdir}/topmed_chr_${chr}.gds
#savenm=${savedir}/topmed_chr_${chr}.out
#savedir=${basedir}/beds/

mkdir -p $savedir


echo Processing $gdsfl
echo Will save to $savenm


# For GDS creation w/ default filters:
#R CMD BATCH --vanilla "--args --gds.file ${gdsfl} --min.AVGDP 10 --filterCat PASS --min.MAF 0.05 --max.miss 0.05 --removeSNPGDS TRUE --prefix.bed ${savenm}" ${basedir}/Seq2BED_wrapper.R ${savedir}/Seq2BED_wrapper_${chr}.Rout

# For regenie, converting full GDS to bed, no filters:
R CMD BATCH --vanilla "--args --gds.file ${gdsfl} --min.AVGDP 0 --filterCat PASS --min.MAF 0 --max.miss 0 --removeSNPGDS TRUE --prefix.bed ${savenm}" ${basedir}/Seq2BED_wrapper.R ${savedir}/Seq2BED_wrapper_${chr}.Rout

