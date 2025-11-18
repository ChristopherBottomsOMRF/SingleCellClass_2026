#!/usr/bin/bash -l
#SBATCH --cpus-per-task 5
#SBATCH --mem 20G

# NOTE: Feel free to copy this file and change it however is helpful to you.

#==============================================================================
# BASH Strict mode (i.e. "fail fast" to reduce hard-to-find bugs)
set -e          # EXIT the script if any command returns non-zero exit status.
set -E          # Make ERR trapping work inside functions too.
set -u          # Variables must be pre-defined before using them.
set -o pipefail # If a pipe fails, returns the error code for the failed pipe
                #  even if it isn't the last command in a series of pipes.

module load cellranger/9.0.0

# If we don't have exactly 4 arguments, give usage message and exit
if [[ "$#" -ne 4 ]]; then
  echo "usage:"
  echo
  echo "  sbatch $0 fastq_id sample_id fastq_dir transcriptome"
  echo
  echo "  fastq_id       base name of FASTQ files (e.g., sample_A for sample_A_S99_L001_R1_001.fastq.gz)"
  echo "  sample_id      'friendly' sample name corresponding to fastq_id"
  echo "  fastq_dir      directory where FASTQ files are found"
  echo "  transcriptome  full path of cellranger transcriptome reference (e.g., /Volumes/shared-refs/cellranger/refdata-gex-GRCh38-2024-A)"
  echo
  echo "  NOTE: This script will not run without using sbatch!"
  exit
fi

# Capture sample ID (that matches the FASTQ filename)
fastq_id=$1
sample_id=$2
fastq_dir=$3
transcriptome=$4

# Get number of CPUS and memory from scheduler variables
# NOTE: This is why you must use sbatch to run this script
cpus=$SLURM_CPUS_PER_TASK
memory=$(( $SLURM_MEM_PER_NODE / 1024 )) # divide SLURM memory (in MB) to get GB
                                         # need 16G+ for use with GRCh38-2020-A (human reference)

# Run cellranger for each sample
#
# Backslashes "\" at the end of each line makes a "line continuation", so the
# following line is effectively just a continuation of the current line.
# Warning: any character (e.g., spaces) after the \ will "break" the line
#
# cellranger --disable-ui prevents the UI interface from being run, which
#
cellranger count \
    --fastqs $fastq_dir \
    --id $sample_id \
    --sample $fastq_id \
    --create-bam true \
    --transcriptome $transcriptome \
    --localcores $cpus \
    --localmem $memory \
    --jobmode local
