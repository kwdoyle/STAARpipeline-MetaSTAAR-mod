#!/bin/bash

#$ -N staarslidwin
#$ -t 1-573
#$ -tc 5
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

scriptdir="/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/"

id=${SGE_TASK_ID}

command=sliding_window

#export R_LIBS_USER=$HOME/R-3.6.1-MKL
#echo $R_LIBS_USER

#/n/home05/zilinli/R-3.6.1/bin/Rscript --slave --no-restore --no-save STAARpipeline_Sliding_Window.r ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout

bash ${scriptdir}/staar_wrapper.sh $command $id

