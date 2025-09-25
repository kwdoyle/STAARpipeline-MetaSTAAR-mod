# STAAR And MetaSTAAR Pipeline

This repo functions as an all-encompasing wrapper around both the STAARpipeline-Tutorial and MetaSTAARlite-Tutorial repositories to run STAAR and MetaSTAAR, respectively.\
The progression of steps run here mirrors the same steps from both of these repositories, easily facilitating them through `staar_wrapper.sh` \
\
Note that most steps of this analysis are meant to be run on a job-submission cluster. The current code is set to run array jobs on a Sun Grid Engine cluster (as seen in `run_staar_wrapper.sh` 
and all scripts in the `batch jobs` directories, but this can be modified to suit whichever cluster of your choosing. \
\
Commands will be submitted to your cluster like so: \
`qsub ./run_staar_wrapper.sh [argument]` \
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Gene_Centric_Noncoding.sh` \
\
Some pre-step commands can be run without job submission: \
`bash ./staar_wrapper.sh prestep` \
\
**Each script within batch_jobs and run_staar_wrapper.sh has a variable pointing to the directory path of this repository, which will need to be updated before running.** \
\
*(Note: there is effectively no difference between run_staar_wrapper.sh and the batch scripts supplied in batch_jobs. \
run_staar_wrapper.sh was created before the scripts from batch_jobs from STAARpipeline-Tutorial and MetaSTAARlite-Tutorial were used)* \
\
Before running `staar_wrapper.sh`, modify the input parameters in this script to correspond to the project name, where to save output to, etc.
Then each of the following steps can be run:

# STAARpipeline Pre-Steps (following those from STAARpipeline-Tutorial):
## Convert vcf (or bed) to gds
`qsub ./run_staar_wrapper.sh gds`

## List variants
`qsub ./run_staar_wrapper.sh listVariants`

## Annotate variants
`qsub ./run_staar_wrapper.sh annotateVariants`

## Annotate gds file (create "agds")
`qsub ./run_staar_wrapper.sh agds`

## Run "prestep"
`bash ./staar_wrapper.sh prestep`

## Fit null model
`bash ./staar_wrapper.sh nullmodel`

# STAARpipeline Analyses (Note: these can be run in any order)
## Individual Variant Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Individual_Analysis.sh`

## Gene-Centric Coding Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Gene_Centric_Coding.sh`

## Gene-Centric Noncoding Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Gene_Centric_Noncoding.sh`

## Gene-Centric ncRNA Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Gene_Centric_ncRNA.sh`

## Sliding Window Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Sliding_Window.sh`

## Dynamic Window Analysis
`qsub ./STAARpipeline/batch_jobs/STAARpipeline_Dynamic_Window.sh`

# MetaSTAAR Analyses 
## (from MetaSTAARlite-Tutorial -- can also be run in any order, although worker functions are first required to run actual meta analysis script)
Each worker script will be run per-cohort, then the main (non-worker) counterpart script will be run across the worker output from each cohort.
## Individual Variant Meta Analysis
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_worker_Individual_Analysis.sh` \
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_Individual_Analysis.sh`

## Gene-Centric Coding Meta Analysis
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_worker_Gene_Centric_Coding.sh` \
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_Gene_Centric_Coding.sh`

## Gene-Centric Noncoding Meta Analysis
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_worker_Gene_Centric_Noncoding.sh` \
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_Gene_Centric_Noncoding.sh`

## Gene-Centric ncRNA Meta Analysis
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_worker_Gene_Centric_ncRNA.sh` \
`qsub ./MetaSTAAR/batch_jobs/MetaSTAARlite_Gene_Centric_ncRNA.sh`

# MetaSTAAR Conditional Analyses
After supplying the list of variants to condition on (known_loci), each conditional script is run via:
## Individual Meta Conditional Analysis
`bash ./staar_wrapper.sh individual_worker_cond`

## Gene-Centric Coding Conditional Analysis
`bash ./staar_wrapper.sh coding_worker_cond` (Note: not yet implemented)

## Gene-Centric Noncoding Conditional Analysis
`bash ./staar_wrapper.sh noncoding_worker_cond`

## Gene-Centric ncRNA Conditional Analysis
`bash ./staar_wrapper.sh ncRNA_worker_cond`
