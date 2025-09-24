#!/bin/bash

#$ -N staardynamwin
#$ -t 1-1879
#$ -tc 7
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

id=${SGE_TASK_ID}

command=dynamic_window

#export R_LIBS_USER=$HOME/R-3.6.1-MKL
#echo $R_LIBS_USER

#/n/home05/zilinli/R-3.6.1/bin/Rscript --slave --no-restore --no-save STAARpipeline_Dynamic_Window.r ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout

bash ~/noncoding_telo/code/shell_scripts/STAAR/run_qsub/qsub_staar_pipeline2.sh $command $id

