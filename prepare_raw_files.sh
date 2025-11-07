#!/bin/bash

#This is a small script to copy the input fastq files as well as the reference genome to the scracth folder of the VSC HPC.
#Meant for reproducibility

#variables
WORKING_DIR="$VSC_SCRATCH/projects/BPK081_subclones"
INPUT_FASTQ="/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_Gabriel_WGS_Ldonovani_20251029/X204SC25056965-Z01-F001.tar"
REFERENCE_GENOME="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/"
BIOINFOTOOLS_DIR="$VSC_SCRATCH/projects/bioinfo_tools/"

#copy fastq files
##go to the project directory
cd "$WORKING_DIR"
#create folder to store raw fastq files
mkdir -p "inputs/WGS/OX_lines"
cd "inputs/WGS/OX_lines"

#copy the data
echo "Coping reference genome..."
mkdir -p inputs/reference_genome
cp -r "$REFERENCE_GENOME" "inputs/reference_genome/"

echo "Coping fastq zipped directory..."
pv "$INPUT_FASTQ" > raw_fastq_files.tar #with `pv` we get a nice progress bar :)

echo "Uncompressing it..."
mkdir raw_fastq_files
pv raw_fastq_files.tar | tar -x -C raw_fastq_files/
mv raw_fastq_files/*/* raw_fastq_files/
rm raw_fastq_files.tar

#check file integrity
echo "Checking file integrity"
cd "raw_fastq_files"
md5sum -c MD5.txt

#move all fastq files to the same directory and delete the rest
echo "Moving all fastq files to the same directory"
find . -type f -iname "*.fq.gz" -exec mv "{}" . \;
find . ! -iname "*.fq.gz" -delete

#rename the files to simpler names
echo "simplifying file names..."
for f in *.fq.gz; do
    prefix=${f%%_*}                   
    last=${f##*_}                    
    new="${prefix}_R${last%.fq.gz}.fastq.gz"
    echo "$f -> $new"
    mv "$f" "$new"
done

#now we use fastp to clean the reads
echo "Cleaning reads with fastp"
cd $WORKING_DIR/

#first we create a tsv file specifying the sample name, read1 and read2 paths
"$BIOINFOTOOLS_DIR/short_read_seq/write_file_list.sh" -i inputs/WGS/OX_lines/raw_fastq_files/ -o temp.tsv

#now we use the tsv file to send fastp as a slurm array
nfiles=$( cat temp.tsv | wc -l )
sbatch --array=1-${nfiles} --wait "$BIOINFOTOOLS_DIR/short_read_seq/run_fastp.sh" \
    -l "temp.tsv" \
    -o "inputs/WGS/OX_lines/cleaned_fastq_files/"

echo Done!

#now create another tsv file specifying the sample name, read1 and read2 paths for the cleaned reads
"$BIOINFOTOOLS_DIR/short_read_seq/"write_file_list.sh -i inputs/WGS/OX_lines/cleaned_fastq_files/ -o temp.tsv
