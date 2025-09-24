# This script is sourced by staar_wrapper.sh when submitted by qsub
# so that the compute nodes know which paths to use

# specify the exact R and r libs I need
alias R="/opt/R/R-4.3.1/bin/R"
#export R_LIBS="/groups/garcia/users/kd2630/R/x86_64-pc-linux-gnu-library/4.3"
export R_LIBS="/efs/garcia/users/kd2630/R/x86_64-pc-linux-gnu-library/4.3:/opt/R/R-4.3.1/lib64/R/library"

# path info from .bash_profile
PATH=$PATH:$HOME/.local/bin:$HOME/bin
# add my own cmake
PATH=~/cmake-3.24.2/bin/:~/libxml2-2.9.9/bin/:$PATH
# adding my own
PATH=/opt/R/R-4.3.1/bin/:$PATH
# adding bcftools
PATH=~/bcftools/:$PATH
# adding htslib (mainly for a tabix that works)
PATH=~/htslib/:$PATH

export PATH

