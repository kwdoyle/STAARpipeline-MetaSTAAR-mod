#!/bin/bash

#$ -N staarcoding
#$ -t 175-179
#$ -tc 5
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

id=${SGE_TASK_ID}

command=coding_analysis

#export R_LIBS_USER=$HOME/R-3.6.1-MKL
#echo $R_LIBS_USER

#/n/home05/zilinli/R-3.6.1/bin/Rscript --slave --no-restore --no-save STAARpipeline_Gene_Centric_Coding.r ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout
#

bash ~/noncoding_telo/code/STAAR/shell_scripts/run_qsub/qsub_staar_pipeline2.sh $command $id


# -t 1-379
