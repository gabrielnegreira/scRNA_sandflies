#!/bin/bash

#This script is used to set common functions and environment variables

#--- Parameters ---
##reference genomes
WGS_REF="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/LdCL_with_kDNA.fasta"
scRNA_REF="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/cellrangerRNA_formatted/"
scDNA_REF="/scratch/antwerpen/grp/aitg/jcdujardin/reference_genomes/L_donovani/LdCL/cellrangerDNA_formatted/"

##working dirs
HPC_WD="${VSC_SCRATCH}/projects/BPK081_subclones/"

##parameters
WGS_BIN_SIZE=1000

#--- functions ---
#create a function to call the `clean_fastq_names.sh` script
clean_fastq_names(){
    echo "[INFO] simplifying the names of fastq files..."
    "${HPC_WD}/git_modules/Illumina_tools/clean_fastq_names.sh" "$@"
}

#create a function to call the `create_file_list.sh` script
file_list(){
    echo "[INFO] listing fastq files..."
    "${HPC_WD}/git_modules/Illumina_tools/write_file_list.sh" "$@"
}
#create a function to run fastp
run_fastp(){
    echo "[INFO] Running fastp on a total of $nfiles files..."
    sbatch --array=1-${nfiles}%20 --wait "${HPC_WD}/git_modules/Illumina_tools/run_fastp.sh" "$@"
    mkdir -p fastp_slurm_logs
    mv "slurm-*" fastp_slurm_logs/
}

#create a function to run BWA-MEM
run_bwa(){
    echo "[INFO] Running BWA-mem to map reads from $nfiles files to the reference genome."
    sbatch --array=1-${nfiles}%20 --wait "${HPC_WD}/git_modules/Illumina_tools/run_BWA-MEM.sh" "$@"
    mkdir -p bwa_slurm_logs
    mv "slurm-*" bwa_slurm_logs/
}

#create a function to calculate binned counts from bam files
binned_counts(){
    echo "[INFO] counting reads for a total of $nfiles bam files."
    sbatch --array=1-${nfiles}%20 --wait "${HPC_WD}/git_modules/WGS_tools/calc_read_count.sh" "$@"
    mkdir -p count_slurm_logs
    mv "slurm-*" count_slurm_logs/
}

#create a function to calculate the mappability of the reference genome
ref_mappability(){
    echo "[INFO] estimating mappability for ${BIN_SIZE}bp bins in the $( basename $REF_GENOME ) reference genome"
    sbatch --wait "${HPC_WD}/git_modules/WGS_tools/compute_gc_and_mappability.sh" "$@"
    mkdir -p mappability_slurm_logs
    mv "slurm-*" mappability_slurm_logs/
}

#create a function to run the aneuploidy analysis
get_aneuploidy(){
    module --force purge
    module load calcua/all calcua/2025a
    module load R/4.5.1-gfbf-2025a
    module load R-bundle-CRAN/2025.10-foss-2025a
    #specifies a path to install R packages on my `data` directory.
    sed -i '/^R_LIBS_USER=/d' "$VSC_HOME/.Renviron"
    echo 'R_LIBS_USER=${VSC_DATA}/Rlibs/${VSC_OS_LOCAL}/${VSC_ARCH_LOCAL}/R-${EBVERSIONR}' >> "$VSC_HOME/.Renviron"
    
    echo "[INFO] determining aneuploidies..."
    Rscript "${HPC_WD}/git_modules/WGS_tools/get_aneuploidy.R" "$@"
}