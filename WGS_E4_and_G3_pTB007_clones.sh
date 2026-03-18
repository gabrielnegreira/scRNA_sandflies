#!/bin/bash

#This script contains the workflow I followed to compute binned counts from fastq files regarding the BPK081 G3 and E4 pTB007 clones.
set -euo pipefail
IFS=$'\n\t'

#source the common variables and functions
source set_env.sh

#set variables
FASTQ_DIR="/data/antwerpen/grp/aitg/jcdujardin/LeishMAP_WGS_E4_subclones_gabriel_20260313/X204SC26024633-Z01-F001/01.RawData/"
OUTPUTS_DIR="$HPC_WD/inputs/WGS/E4_G3_pTB007/"

#run the WGS pipeline
rm -rf "$OUTPUTS_DIR"
./WGS_pipeline.sh \
    --fastq-dir "$FASTQ_DIR" \
    --outputs-dir "$OUTPUTS_DIR" \
    --ref-genome "$WGS_REF" \
    --bin-size 1000

