#!/bin/bash
#This is the full pipeline used to performe the WGS analysis. Basically from fastq to aneuploidy determination.
# --- preflight --- 
##set usage helper
set -euo pipefail

usage() {
  cat <<EOF
Usage: sbatch $(basename "$0") [options]

Required:
  --fastq-dir   <path>   Path to the fastq files.
  --outputs-dir <path>   Path to the directory where the outputs will be stored.
  --ref-genome  <path>   Path to the reference genome fasta file.


Optional:
  --bin-size    <int>    Bin size in bp for read counting. Defaults to 20000.
  -h, --help             Show this help
EOF
  exit 1
}

# --- arguments ---

# Defaults
FASTQ_DIR=""
OUTPUTS_DIR=""
REF_GENOME=""
BIN_SIZE="20000"

# Parse
while [[ $# -gt 0 ]]; do
  case $1 in
    --fastq-dir)
      FASTQ_DIR="$2"
      shift 2
      ;;
    --outputs-dir)
      OUTPUTS_DIR="$2"
      shift 2
      ;;
    --ref-genome)
      REF_GENOME="$2"
      shift 2
      ;;
    --bin-size)
      BIN_SIZE="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "ERROR: Unknown option '$1'" >&2
      echo "Use --help for usage information" >&2
      exit 1
      ;;
  esac
done

# Validate
if [[ -z "$FASTQ_DIR" || -z "$OUTPUTS_DIR" || -z "$REF_GENOME" ]]; then
  echo "ERROR: --fastq-dir, --output-dir, --ref-genome and --bin-size are required" >&2
  usage
fi

echo ">>> Parameters:"
echo "    Fastq directory:   $FASTQ_DIR"
echo "    Output dir:        $OUTPUTS_DIR"
echo "    Reference genome:  $REF_GENOME"
echo "    Bin size:          ${BIN_SIZE}bp"
echo


# --- actual code ---
#source the common variables and functions
script_path="$(realpath "$(dirname "$0")")"
source "${script_path}/set_env.sh"

#create folder to store raw fastq files
mkdir -p "${OUTPUTS_DIR}/fastq_files"
cd "${OUTPUTS_DIR}"

#copy the fastq files
echo "[INFO] Coping fastq files"
find "$FASTQ_DIR" -type f -regextype posix-extended -regex '.*\.(fastq|fq)(\.gz)?$' -exec cp -u {} "fastq_files/" \;

#change permissions
chmod -R 770 "fastq_files/"

#rename the files to simpler names
clean_fastq_names --fastq_dir "fastq_files/"

#now that file names are clean, we write a tsv file listing them to be able to submit them as slurm arrays.
file_list -i "fastq_files/" -o temp.tsv

#now we use the tsv file to send fastp as a slurm array
nfiles=$( cat temp.tsv | wc -l )
run_fastp \
    -l "temp.tsv" \
    -o "cleaned_fastq_files/"

#now create another tsv file specifying the sample name, read1 and read2 paths for the cleaned reads
rm temp.tsv
file_list -i "cleaned_fastq_files/" -o temp.tsv

#use it to run BWA-MEM 
nfiles=$( cat temp.tsv | wc -l )
run_bwa \
    -l temp.tsv \
    -r $REF_GENOME \
    -o "bam_files/"

#now use the generated bam files to count binned reads
# Should be:
find "bam_files/" -maxdepth 1 -name "*.bam" > temp.tsv

nfiles=$( cat temp.tsv | wc -l )
binned_counts \
    -l temp.tsv \
    -r "$REF_GENOME" \
    -o "count_files/" \
    -s $BIN_SIZE

#finally we estimate genome mappability and bind it to the output count files
rm -f *mappability.temp
ref_mappability \
    --ref "$REF_GENOME" \
    --output-dir . \
    --bin-size $BIN_SIZE

map_file=$( echo *mappability.temp )

#then we merge the files
tail -n +2 "$map_file" | cut -f9 > mappability_column.temp

for bed_file in "${OUTPUTS_DIR}/count_files/"*.bed; do
    paste "$bed_file" mappability_column.temp > "${bed_file}.temp"
    mv "${bed_file}.temp" "$bed_file"
done

#remove temporary files
rm *.temp

#now perform aneuploidy analysis
get_aneuploidy \
    --inputs_dir="count_files/" \
    --outputs_dir="aneuploidy_analysis/" \
    --ploidy=2

echo "[DONE] Files stored in ${OUTPUTS_DIR}/"