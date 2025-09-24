#!/bin/bash
source /efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/makeenv2.sh

# This script sources my env and actually runs the conversion script

chr=$1
gdsdir=~/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/
basedir=~/noncoding_telo/STAAR/TOPMed_Full_Cohort_grm/GRM/
savedir=${basedir}/beds/

mkdir -p $savedir

gdsfl=${gdsdir}/topmed_chr_${chr}.gds
savenm=${savedir}/topmed_chr_${chr}.out

echo Processing $gdsfl
echo Will save to $savenm


R CMD BATCH --vanilla "--args --gds.file ${gdsfl} --min.AVGDP 10 --filterCat PASS --min.MAF 0.05 --max.miss 0.05 --removeSNPGDS TRUE --prefix.bed ${savenm}" ${basedir}/Seq2BED_wrapper.R ${savedir}/Seq2BED_wrapper_${chr}.Rout

