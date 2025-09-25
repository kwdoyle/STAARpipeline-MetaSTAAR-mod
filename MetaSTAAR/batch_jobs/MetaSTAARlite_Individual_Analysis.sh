#!/bin/bash

#$ -N metastaarindiv
#$ -t 1-294
#$ -tc 7
#$ -V
#$ -cwd
#$ -l h_vmem=20G
#$ -j yes
#$ -q all.q

id=${SGE_TASK_ID}

command=individual_meta

#module purge
#module load r/4.3.1
#export R_LIBS_USER=/proj/xihaoli_lab/users/xihaoli/R-4.3.1
#R --vanilla --args ${SLURM_ARRAY_TASK_ID} < $1 > "${1}.${SLURM_ARRAY_TASK_ID}.out" 2> "${1}.${SLURM_ARRAY_TASK_ID}.err"

bash /efs/garcia/users/kd2630/noncoding_telo/code/STAAR/shell_scripts/run_qsub/qsub_staar_pipeline2_metastaar.sh $command $id

# TODO not all of these ran for some reason: 135,233,238,275,291
# oookkayy, it's because of some deep core issue with this indiv analysis pipeline where I see this error:
# longer object length is not a multiple of shorter object length
# $ -t 1-294

