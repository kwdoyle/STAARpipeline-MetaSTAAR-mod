#!/bin/bash

# use SGE directives
#$ -t 1-22
#$ -tc 5
#$ -V
#$ -l mem_free=50G
#$ -l s_vmem=50G
#$ -l h_vmem=50G
#$ -b yes
#$ -cwd
#$ -j yes
#$ -q all.q
#$ -N staar

chr=${SGE_TASK_ID}

# NOTE: this script is used to run all the pre-analysis steps in the staar pipeline that require job ids

# task to run
command=agds  #agds  #annotateVariants  #listVariants  #indiv_known_loci  #individual_analysis  #coding_analysis  #ncRNA_analysis  #noncoding_analysis  #gds

bash ~/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/staar_wrapper.sh $command $chr

