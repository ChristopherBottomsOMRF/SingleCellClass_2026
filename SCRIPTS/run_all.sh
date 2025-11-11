#!/usr/bin/bash -l 
#SBATCH --cpus-per-task 30
#SBATCH --mem-per-cpu 4G

if [ "$SLURM_JOB_ID" == "" ]; then
    # Capture name of the directory this script is located in and add it to PATH
    export script_dir=$(dirname $(realpath $0))
    export PATH="$script_dir:$PATH"
    #------------------------------------------------------------------------------

    exec sbatch $0 $@ 
fi

#==============================================================================
# BASH Strict mode (i.e. "fail fast" to reduce hard-to-find bugs)
set -e          # EXIT the script if any command returns non-zero exit status.
set -E          # Make ERR trapping work inside functions too.
set -u          # Variables must be pre-defined before using them.
set -o pipefail # If a pipe fails, returns the error code for the failed pipe
                #  even if it isn't the last command in a series of pipes.

module load R/4.5.1-mkl seurat/5.3

# create_prefiltered_merged_seurat_obj.R
# Script_for_selecting_filter_thresholds.server_version.R
# StandardProcessing_DoubletRemoval_and_Integration.R
Plot_integrations.R
