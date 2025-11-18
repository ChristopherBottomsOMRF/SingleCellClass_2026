#!/usr/bin/bash -l

#===============================================================================
# Allow running other scripts within the same directory as this script
#==========

# Capture name of the directory this script is located in and add it to PATH
export script_dir=$(dirname $(realpath $0))
export PATH="$script_dir:$PATH"
#------------------------------------------------------------------------------

# If we don't have exactly 3 arguments, give usage message and exit
if [[ "$#" -ne 3 ]]; then
  echo "usage:"
  echo
  echo "  $0 sample_sheet fastq_dir transcriptome"
  echo
  echo "  sample_sheet   CSV file with a header, containing 3 columns: Fastq_ID,group,replicate"
  echo "  fastq_dir      directory where FASTQ files are found"
  echo "  transcriptome  full path of cellranger transcriptome reference (e.g., /Volumes/shared-refs/cellranger/refdata-gex-GRCh38-2024-A)"
  exit
fi

sample_sheet=$1

# Export variables so that parallel can use them.
export fastq_dir=$2
export transcriptome=$3

# Run cellranger

module load slurm

cat $sample_sheet |
  tail -n +2 |      # skip the header line
  parallel --colsep=, 'sbatch --job-name {2}_{3}_cellranger cellranger_wrapper.sh {1} {2}_{3} $fastq_dir $transcriptome'
