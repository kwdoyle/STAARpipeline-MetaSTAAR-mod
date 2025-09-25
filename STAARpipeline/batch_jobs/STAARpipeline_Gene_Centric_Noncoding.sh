#!/bin/bash

#$ -N staarnoncoding
#$ -t 1-379
#$ -tc 5
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

scriptdir="/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/"

id=${SGE_TASK_ID}

command=noncoding_analysis

#export R_LIBS_USER=$HOME/R-3.6.1-MKL
#echo $R_LIBS_USER

#/n/home05/zilinli/R-3.6.1/bin/Rscript --slave --no-restore --no-save STAARpipeline_Gene_Centric_Noncoding.r ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout

#bash ~/noncoding_telo/code/shell_scripts/STAAR/run_qsub/qsub_staar_pipeline2.sh $command $id
bash ${scriptdir}/staar_wrapper.sh $command $id

