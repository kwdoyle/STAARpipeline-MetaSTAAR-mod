#!/bin/bash

# use SGE directives
#$ -N seq2bed
#$ -t 1-22
#$ -tc 5
#$ -V
#$ -b yes
#$ -j yes
#$ -cwd
#$ -l h_vmem=10G
#$ -q all.q

scriptdir="/efs/garcia/users/kd2630/noncoding_telo/STAAR/STAARpipeline-MetaSTAAR-mod/"

chr=${SGE_TASK_ID}

bash ${scriptdir}/Seq2BED_wrapper.sh ${chr}

