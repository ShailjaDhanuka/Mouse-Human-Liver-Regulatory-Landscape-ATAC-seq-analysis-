#!/usr/bin/env bash
#SBATCH --job-name=open_chromatin
#SBATCH --partition=RM
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Usage: sbatch open_chromatin_finding_script.sh \
#          <MouseToHuman.HALPER.narrowPeak.gz> \
#          <HumanToMouse.HALPER.narrowPeak.gz> \
#          <mouse_atac.narrowPeak> \
#          <human_atac.narrowPeak> \
#          <output_dir>

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <MouseToHuman.HALPER.narrowPeak.gz> <HumanToMouse.HALPER.narrowPeak.gz> <mouse_atac.narrowPeak> <human_atac.narrowPeak> <output_dir>"
    exit 1
fi

MOUSE_TO_HUMAN_GZ="$1"
HUMAN_TO_MOUSE_GZ="$2"
MOUSE_ATAC="$3"
HUMAN_ATAC="$4"
OUTDIR="$5"

mkdir -p "$OUTDIR"

module load bedtools

# 1. prepare BED
zcat "$MOUSE_TO_HUMAN_GZ" | sort -k1,1 -k2,2n > "$OUTDIR/mouse_to_human_orthologs.bed"
zcat "$HUMAN_TO_MOUSE_GZ" | sort -k1,1 -k2,2n > "$OUTDIR/human_to_mouse_orthologs.bed"
zcat "$MOUSE_ATAC" | sort -k1,1 -k2,2n > "$OUTDIR/mouse_atac.bed"
zcat "$HUMAN_ATAC" | sort -k1,1 -k2,2n > "$OUTDIR/human_atac.bed"

# 2. shared open
bedtools intersect -a "$OUTDIR/mouse_to_human_orthologs.bed" -b "$OUTDIR/human_atac.bed" -u -sorted > "$OUTDIR/mouse_shared_in_human.humanCoords.bed"
bedtools intersect -a "$OUTDIR/human_to_mouse_orthologs.bed" -b "$OUTDIR/mouse_atac.bed" -u -sorted > "$OUTDIR/human_shared_in_mouse.mouseCoords.bed"

# 3. species-specific but ortholog exists
bedtools intersect -a "$OUTDIR/mouse_to_human_orthologs.bed" -b "$OUTDIR/human_atac.bed" -v -sorted > "$OUTDIR/mouse_open_human_closed.humanCoords.bed"
bedtools intersect -a "$OUTDIR/human_to_mouse_orthologs.bed" -b "$OUTDIR/mouse_atac.bed" -v -sorted > "$OUTDIR/human_open_mouse_closed.mouseCoords.bed"

# 4. unmappable / no ortholog
comm -23 <(cut -f4 "$OUTDIR/mouse_atac.bed" | sort) \
         <(cut -f4 "$OUTDIR/mouse_to_human_orthologs.bed" | sort) | \
grep -Fwf - "$OUTDIR/mouse_atac.bed" > "$OUTDIR/mouse_no_ortholog.mouseCoords.bed"

comm -23 <(cut -f4 "$OUTDIR/human_atac.bed" | sort) \
         <(cut -f4 "$OUTDIR/human_to_mouse_orthologs.bed" | sort) | \
grep -Fwf - "$OUTDIR/human_atac.bed" > "$OUTDIR/human_no_ortholog.humanCoords.bed"