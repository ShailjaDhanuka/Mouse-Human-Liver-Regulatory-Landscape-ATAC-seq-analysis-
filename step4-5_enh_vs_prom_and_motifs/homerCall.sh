#!/bin/bash

# extra wrapper script just for calling main wrapper script without having to reput all the parameters during testing

#SBATCH --job-name=homer_motif
#SBATCH --partition=RM-shared
#SBATCH --mem=8000MB
#SBATCH --cpus-per-task=4
#SBATCH --time=07:00:00

#Usage: enhancer_vs_promoter.sh \
#   <mouse_peaks_shared_human> \
#   <human_peaks_shared_mouse> \
#   <unique_mouse_human_closed> \
#   <unique_human_mouse_closed> \
#   <mouse_only_no_ortho> \
#   <human_only_no_ortho> \
#   <output_dir>

module load anaconda3
conda activate homer_env
module load bedtools

MOUSE_SHARED_HUMAN="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/mouse_shared_in_human.humanCoords.bed"
HUMAN_SHARED_MOUSE="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/human_shared_in_mouse.mouseCoords.bed"
MOUSE_HUMAN_CLOSED="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/mouse_open_human_closed.humanCoords.bed"
HUMAN_MOUSE_CLOSED="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/human_open_mouse_closed.mouseCoords.bed"
MOUSE_NO_ORTHO="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/mouse_no_ortholog.mouseCoords.bed"
HUMAN_NO_ORTHO="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result/human_no_ortholog.humanCoords.bed"
OUTPUT_DIR="/ocean/projects/bio230007p/sturecki/project/HOMER_run"

bash enhancer_vs_promoter.sh \
    $MOUSE_SHARED_HUMAN \
    $HUMAN_SHARED_MOUSE \
    $MOUSE_HUMAN_CLOSED \
    $HUMAN_MOUSE_CLOSED \
    $MOUSE_NO_ORTHO \
    $HUMAN_NO_ORTHO \
    $OUTPUT_DIR

