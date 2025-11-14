#!/bin/bash

set -euo pipefail
IFS=$'\n\t' #this sets shell split strings to only tabs and newlines. 

#variables
WORKING_DIR="$VSC_SCRATCH/projects/BPK081_subclones"
INPUT_FASTQ="/data/antwerpen/grp/aitg/jcdujardin/MADLEI_wgs_20210610/"
REFERENCE_GENOME="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/"
BIOINFOTOOLS_DIR="$VSC_SCRATCH/projects/bioinfo_tools/"

#copy fastq files
##go to the project directory
cd "$WORKING_DIR"

#copy the reference genome
echo "Coping reference genome..."
mkdir -p inputs/reference_genome
cp -r "$REFERENCE_GENOME" "inputs/reference_genome/"

#copy input fastq
echo copying fastq files to scratch...
mkdir -p inputs/WGS/subclones/raw_fastq_files
find "$INPUT_FASTQ" \
  -type f \
  -regextype posix-extended \
  -regex '.*104487-001-(04[1-9]|0[5-9][0-9]|10[0-3]).*' \
  -exec cp {} "inputs/WGS/subclones/raw_fastq_files/" \;

#in this specific dataset, read2 is actually named as `R3`, as R2 is for indexes. So we remove R2 files and rename R3 into R2.
echo deleting index files...
cd inputs/WGS/subclones/raw_fastq_files
rm -f *R2.fastq.gz
echo renaming R3 files into R2...
for f in *.fastq.gz; do
    if [[ "$f" == *R3.fastq.gz ]]; then
        mv -v -- "$f" "${f/R3.fastq.gz/R2.fastq.gz}"
    fi    
done

#change file names to make it usable with the `write_file_list.sh` script.
echo "replacing first underscore in file names..."
for f in *.fastq.gz; do
  [[ "$f" == *_* ]] || continue
  mv -v -- "$f" "${f/_/-}"
done

#rename the files to simpler names
echo "simplifying file names..."
for f in *.fastq.gz; do
    prefix=${f%%_*}                   
    last=${f##*_}                    
    new="${prefix}_${last%.fastq.gz}.fastq.gz"
    echo "$f -> $new"
    mv "$f" "$new"
done

#create a temporary list file to be used with fastp later
cd "$WORKING_DIR"
"$BIOINFOTOOLS_DIR/short_read_seq/write_file_list.sh" -i inputs/WGS/subclones/raw_fastq_files/ -o temp.tsv

#use it to run fastp as a slurm array
#now we use the tsv file to send fastp as a slurm array
nfiles=$( cat temp.tsv | wc -l )
echo "running fastp on a total of $nfiles files."
sbatch --array=1-${nfiles}%16 --wait "$BIOINFOTOOLS_DIR/short_read_seq/run_fastp.sh" \
    -l "temp.tsv" \
    -o "inputs/WGS/subclones/cleaned_fastq_files/"

#create another temporary list file to be used with BWA MEM now.
"$BIOINFOTOOLS_DIR/short_read_seq/write_file_list.sh" -i inputs/WGS/subclones/cleaned_fastq_files/ -o temp.tsv

#use it to run BWA-MEM 
nfiles=$( cat temp.tsv | wc -l )
echo "running bwa mem on a total of $nfiles files..."
sbatch --array=1-${nfiles}%14 --wait "$BIOINFOTOOLS_DIR/short_read_seq/run_BWA-MEM.sh" \
    -l temp.tsv \
    -r "inputs/reference_genome/LdCL_original.fasta" \
    -o "inputs/WGS/subclones/bam_files/"

#now use the generated bam files to count binned reads
find "inputs/WGS/subclones/bam_files/"*.bam > temp.tsv 
nfiles=$( cat temp.tsv | wc -l )
echo "counting binned reads for a total of $nfiles files..."
sbatch --array=1-${nfiles}%14 --wait "$BIOINFOTOOLS_DIR/short_read_seq/calc_read_count.sh" \
    -l temp.tsv \
    -r "inputs/reference_genome/LdCL_original.fasta" \
    -o "inputs/WGS/subclones/count_files/" \
    -s 1000
