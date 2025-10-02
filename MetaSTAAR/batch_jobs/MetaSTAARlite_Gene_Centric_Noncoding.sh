#!/bin/bash

#$ -N metastaarnoncoding
#$ -t 1-379
#$ -tc 7
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

scriptdir="/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/"

id=${SGE_TASK_ID}

command=noncoding_meta

#module purge
#module load r/4.3.1
#export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.3.1
#R --vanilla --args ${SLURM_ARRAY_TASK_ID} < $1 > "${1}.${SLURM_ARRAY_TASK_ID}.out" 2> "${1}.${SLURM_ARRAY_TASK_ID}.err"

bash ${scriptdir}/staar_wrapper.sh $command $id

