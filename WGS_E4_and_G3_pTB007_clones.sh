#!/bin/bash

#This script contains the workflow I followed to compute binned counts from fastq files regarding the BPK081 G3 and E4 pTB007 clones.

set -euo pipefail
IFS=$'\n\t'

#variables
WORKING_DIR="${VSC_SCRATCH}/projects/BPK081_subclones/"
INPUT_FASTQ="/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_WGS_E4_subclones_gabriel_20260313/X204SC26024633-Z01-F001/01.RawData/"
REF_GENOME="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/LdCL_original.fasta"
BIN_SIZE=1000

#functions
#create a function to call the `create_file_list.sh` script
file_list(){
    echo "[INFO] listing fastq files..."
    "${WORKING_DIR}/git_modules/Illumina_tools/write_file_list.sh" "$@"
}
#create a function to run fastp
run_fastp(){
    echo "[INFO] Running fastp on a total of $nfiles files..."
    sbatch --array=1-${nfiles}%20 --wait "${WORKING_DIR}/git_modules/Illumina_tools/run_fastp.sh" "$@"
}

#create a function to run BWA-MEM
run_bwa(){
    echo "[INFO] Running BWA-mem to map reads from $nfiles files to the reference genome."
    sbatch --array=1-${nfiles}%20 --wait "${WORKING_DIR}/git_modules/Illumina_tools/run_BWA-MEM.sh" "$@"
}

#create a function to calculate binned counts from bam files
binned_counts(){
    echo "[INFO] counting reads for a total of $nfiles bam files."
    sbatch --array=1-${nfiles}%20 --wait "${WORKING_DIR}/git_modules/WGS_tools/calc_read_count.sh" "$@"
}

#create a function to calculate the mappability of the reference genome
ref_mappability(){
    echo "[INFO] estimating mappability for ${BIN_SIZE}bp bins in the $( basename $REF_GENOME ) reference genome"
    sbatch --wait "${WORKING_DIR}/git_modules/WGS_tools/compute_gc_and_mappability.sh" "$@"
}

#create a function to run the aneuploidy analysis
get_aneuploidy(){
    module --force purge
    module load calcua/all calcua/2025a
    module load R/4.5.1-gfbf-2025a
    module load R-bundle-CRAN/2025.10-foss-2025a
    echo "[INFO] determining aneuploidies..."
    Rscript "${WORKING_DIR}/git_modules/WGS_tools/get_aneuploidy.R" "$@"
}

#copy fastq files
##go to the project directory
cd "$WORKING_DIR"

#copy the fastq files
echo "[INFO] Coping fastq files"
#create folder to store raw fastq files
mkdir -p "inputs/WGS/E4_G3_pTB007"
cd "inputs/WGS/E4_G3_pTB007"
find "$INPUT_FASTQ" -name *.fq.gz -print0 -exec cp -u {} . \;

#change permissions
chmod -R 700 .
#rename the files to simpler names
echo "[INFO] simplifying file names..."
for f in *.fq.gz; do
  echo $f
   prefix=${f%%_*}                   
   last=${f##*_}                    
   new="${prefix}_R${last%.fq.gz}.fastq.gz"
   echo "$f -> $new"
   mv "$f" "$new"
done

#now we use fastp to clean the reads
echo "[INFO] Cleaning reads with fastp"
cd "$WORKING_DIR"

 #first we create a tsv file specifying the sample name, read1 and read2 paths for the raw reads
file_list -i inputs/WGS/E4_G3_pTB007/ -o temp.tsv

#now we use the tsv file to send fastp as a slurm array
nfiles=$( cat temp.tsv | wc -l )
run_fastp \
    -l "temp.tsv" \
    -o "inputs/WGS/E4_G3_pTB007/cleaned_fastq_files/"

#now create another tsv file specifying the sample name, read1 and read2 paths for the cleaned reads
rm temp.tsv
file_list -i inputs/WGS/E4_G3_pTB007/cleaned_fastq_files/ -o temp.tsv

#use it to run BWA-MEM 
nfiles=$( cat temp.tsv | wc -l )
run_bwa \
    -l temp.tsv \
    -r $REF_GENOME \
    -o "inputs/WGS/E4_G3_pTB007/bam_files/"

#now use the generated bam files to count binned reads
find "inputs/WGS/E4_G3_pTB007/bam_files/"*.bam > temp.tsv 
nfiles=$( cat temp.tsv | wc -l )
binned_counts \
    -l temp.tsv \
    -r "inputs/reference_genome/LdCL_original.fasta" \
    -o "inputs/WGS/E4_G3_pTB007/count_files/" \
    -s $BIN_SIZE

#finally we estimate genome mappability and bind it to the output count files
rm -f *mappability.tsv
ref_mappability \
    --ref "$REF_GENOME" \
    --output-dir . \
    --bin-size $BIN_SIZE

map_file=$( echo *mappability.tsv )

#then we merge the files
tail -n +2 "$map_file" | cut -f9 > mappability_column.temp

for bed_file in "inputs/WGS/E4_G3_pTB007/count_files/"*.bed; do
    paste "$bed_file" mappability_column.temp > "${bed_file}.tmp"
    mv "${bed_file}.tmp" "$bed_file"
done

#remove temporary files
rm *.temp

echo "[DONE] Files stored in inputs/WGS/E4_G3_pTB007/"

#now perform aneuploidy analysis
get_aneuploidy \
    --inputs_dir="inputs/WGS/E4_G3_pTB007/count_files/" \
    --outputs_dir="inputs/WGS/E4_G3_pTB007/aneuploidy_analysis/" \
    --ploidy=2