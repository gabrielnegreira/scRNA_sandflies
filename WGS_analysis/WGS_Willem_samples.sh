#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=WGS_HAP2
#SBATCH --mail-type=BEGIN,END,FAIL

#This script contains the workflow I followed to compute binned counts from fastq files regarding the HAP2 OX and HAP2.2 KO lines from Willem.

#----- Set the environment ------
set -euo pipefail
IFS=$'\n\t'
source set_env.sh

#----- set variables -----
FASTQ_TAR_FILE="/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_Gabriel_WGS_Willem_samples_20260611/X205SC26056316-Z01-F001.tar" #tar file containing all fastq files downloaded from Novogene
FASTQ_WD="$HPC_WD/inputs/WGS/HAP2_lines" #directory where the fastq files will be stored and processed

#----- Format fastq files -----

#for the pipeline to work we need to first format the fastq file to match the following criteria:
## 1) read 1 files should contain '_R1' in their file name.
## 2) read 2 files should contain '_R2' in their file name.
## 3) Each sample should have only one read 1 and one read 2 file. Merge libraries sequenced in multiple lanes/sequencers!
## 4) file extension should be 'fastq.gz', not 'fq.gz'
## 5) sample name should be specified in the first part of the file name, preceded by an underscore ('_') and not containing an underscore. 
##      example: sample1_R1.fastq
##               sample2_R2.fastq        

wd="$(pwd)" #save current working directory to come back later
startover=false
if $startover; then
    #copy the fastq files
    echo "Coping fastq zipped directory..."
    #create folder to store raw fastq files
    rm -rf "$FASTQ_WD"
    mkdir -p "$FASTQ_WD"
    cd "$FASTQ_WD"
    pv "$FASTQ_TAR_FILE" > "$FASTQ_WD"/raw_fastq_files.tar #with `pv` we get a nice progress bar :)

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

    #fix names
    ## her I am making sure the sample names do not contain underscores.
    echo "fixing file names..."
    for f in *.fq.gz; do
        prefix="${f%%MKDN*}"           # part before MKDN
        rest="${f#*$prefix}"           # from MKDN onward
        new="${prefix//_/}_${rest}"     # remove all _ from prefix, keep rest unchanged
        echo "$f -> $new"
        [ "$f" != "$new" ] && mv "$f" "$new"
    done

    ## now remove anything between the sample name and the lane indicator (since these libraries were run in multiple lanes)
    echo "simplifying file name..."
    for f in *.fq.gz; do
        base="${f%%_L*}"                    # part before L
        prefix="${base%%_*}"                # part before first _
        rest="${f#${base}}"
        new="${prefix}${rest}"
        echo "$f -> $new"
        [ "$f" != "$new" ] && mv "$f" "$new"
    done


    ## Willem mistakenly named his KO lines as OX. Fix it here:
    echo "Fixing KO names..."
    for f in HAP22OX*.fq.gz; do
        new="${f/HAP22OX/HAP2-2KO}"
        echo "$f -> $new"
        mv "$f" "$new"
    done

    ## finally, merge the files of each sample
    echo "merging files from multiple lanes into a single file"
    for sample in $(ls *.fq.gz | sed -E 's/_L[0-9]+_[12]\.fq\.gz$//' | sort -u); do
        r1_files=$(ls ${sample}_L*_1.fq.gz | sort -V)
        r2_files=$(ls ${sample}_L*_2.fq.gz | sort -V)

        echo "Sample: $sample"
        echo "  R1 -> ${sample}_R1.fastq.gz"
        printf '    %s\n' $r1_files
        echo "  R2 -> ${sample}_R2.fastq.gz"
        printf '    %s\n' $r2_files
        echo

        cat $r1_files > ${sample}_R1.fastq.gz &&
        cat $r2_files > ${sample}_R2.fastq.gz &&
        rm $r1_files &&
        rm $r2_files
    done

    ## make sure read pair match was kept
    echo "Checking read pair mates..."
    for r1 in *_R1.fastq.gz; do
        r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"

        if [ ! -f "$r2" ]; then
            echo "$r1 -> missing mate file"
            continue
        fi

        echo "Checking $r1 <-> $r2"

        paste \
        <(zcat "$r1" | sed -n '1~4p' | sed 's/\/1$//' | sed 's/ 1:N:.*$//') \
        <(zcat "$r2" | sed -n '1~4p' | sed 's/\/2$//' | sed 's/ 2:N:.*$//') \
        | awk -v s="$r1" '
            $1 != $2 {print "[FAIL] " s " first mismatch:", $1, $2; bad=1; exit}
            END {
                if (!bad) print "[OK]   " s
            }'
    done
fi

#----- Run the WGS pipeline ---------
OUTPUTS_DIR="$FASTQ_WD/outputs/"
rm -rf "$OUTPUTS_DIR"
cd "$wd"
./WGS_pipeline.sh \
    --fastq-dir "$FASTQ_WD/raw_fastq_files" \
    --outputs-dir "$OUTPUTS_DIR" \
    --ref-genome "$WGS_REF" \
    --bin-size 1000

