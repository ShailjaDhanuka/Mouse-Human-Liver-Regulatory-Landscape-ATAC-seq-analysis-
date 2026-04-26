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

# fail safe - checks if module exists as command in subshell made by script -- skips if it doesn't to prevent error
command -v module &>/dev/null && module load bedtools

# 1. prepare BED
zcat "$MOUSE_TO_HUMAN_GZ" | sort -k1,1 -k2,2n > "$OUTDIR/mouse_to_human_orthologs.bed"
zcat "$HUMAN_TO_MOUSE_GZ" | sort -k1,1 -k2,2n > "$OUTDIR/human_to_mouse_orthologs.bed"
zcat "$MOUSE_ATAC" | sort -k1,1 -k2,2n > "$OUTDIR/mouse_atac.bed"
zcat "$HUMAN_ATAC" | sort -k1,1 -k2,2n > "$OUTDIR/human_atac.bed"

# 2. shared open
bedtools intersect -a "$OUTDIR/mouse_to_human_orthologs.bed" -b "$OUTDIR/human_atac.bed" -u -sorted > "$OUTDIR/mouse_shared_in_human.humanCoords.bed"
bedtools intersect -a "$OUTDIR/human_to_mouse_orthologs.bed" -b "$OUTDIR/mouse_atac.bed" -u -sorted > "$OUTDIR/human_shared_in_mouse.mouseCoords.bed"

# 3. species-specific but ortholog exists (and is closed)

### mouse open, human closed -- keeping mouse coordinates

## finding mouse ATAC regions whose liftover does not overlap human ATAC
bedtools intersect \
  -a "$OUTDIR/mouse_to_human_orthologs.bed" \
  -b "$OUTDIR/human_atac.bed" \
  -v -sorted > "$OUTDIR/human_closed_orthologs.bed"

## getting mouse coordinates whose human orthologs are closed
# HALPER has source (mouse) coords in col 4 as chrom:start-end:summit when input name is "."
# parsing coords, splitting by : and - to get temporary 3 col bed file of mouse coords and intersect with mouse ATAC
awk '{split($4,a,":"); split(a[2],b,"-"); if(length(b)>=2) print a[1]"\t"b[1]"\t"b[2]}' \
  "$OUTDIR/human_closed_orthologs.bed" | sort -k1,1 -k2,2n > "$OUTDIR/tmp_mouse_closed_coords.bed"

## mapping back to original mouse coordinates to get species specific OCRs with closed ortholog
bedtools intersect \
  -a "$OUTDIR/mouse_atac.bed" \
  -b "$OUTDIR/tmp_mouse_closed_coords.bed" \
  -u -sorted > "$OUTDIR/mouse_open_human_closed.mouseCoords.bed"

### human open, mouse closed -- keeping human coordinates

## finding human ATAC regions whose liftover does not overlap mouse ATAC
bedtools intersect \
  -a "$OUTDIR/human_to_mouse_orthologs.bed" \
  -b "$OUTDIR/mouse_atac.bed" \
  -v -sorted > "$OUTDIR/mouse_closed_orthologs.bed"

## getting human coordinates whose mouse orthologs are closed
# HALPER encodes source (human) coords in col 4 as chrom:start-end:summit
# parsing coords, splitting by : and - to get temporary 3 col bed file of human coords, and intersecting with human ATAC
awk '{split($4,a,":"); split(a[2],b,"-"); if(length(b)>=2) print a[1]"\t"b[1]"\t"b[2]}' \
  "$OUTDIR/mouse_closed_orthologs.bed" | sort -k1,1 -k2,2n > "$OUTDIR/tmp_human_closed_coords.bed"

## mapping back to original human coordinates to get species specific OCRs with closed ortholog
bedtools intersect \
  -a "$OUTDIR/human_atac.bed" \
  -b "$OUTDIR/tmp_human_closed_coords.bed" \
  -u -sorted > "$OUTDIR/human_open_mouse_closed.humanCoords.bed"


# 4. unmappable / no ortholog
# peaks whose coordinates do not appear as a source region in any HALPER output row
# extract source coords from col 4 of HALPER output, then use bedtools -v to find ATAC peaks with no orthologs

## parsing coords, splitting by : and - to get temporary 3 col bed file of mouse coords with orthologs
awk '{split($4,a,":"); split(a[2],b,"-"); if(length(b)>=2) print a[1]"\t"b[1]"\t"b[2]}' \
  "$OUTDIR/mouse_to_human_orthologs.bed" | sort -k1,1 -k2,2n > "$OUTDIR/tmp_mouse_has_ortholog_coords.bed"

## negating intersection of mouse coords with orthologs and mouse atac to get unmapped atac peaks
bedtools intersect \
  -a "$OUTDIR/mouse_atac.bed" \
  -b "$OUTDIR/tmp_mouse_has_ortholog_coords.bed" \
  -v -sorted > "$OUTDIR/mouse_no_ortholog.mouseCoords.bed"

## parsing coords, splitting by : and - to get temporary 3 col bed file of human coords with orthologs
awk '{split($4,a,":"); split(a[2],b,"-"); if(length(b)>=2) print a[1]"\t"b[1]"\t"b[2]}' \
  "$OUTDIR/human_to_mouse_orthologs.bed" | sort -k1,1 -k2,2n > "$OUTDIR/tmp_human_has_ortholog_coords.bed"

## negating intersection of human coords with orthologs and human atac to get unmapped atac peaks
bedtools intersect \
  -a "$OUTDIR/human_atac.bed" \
  -b "$OUTDIR/tmp_human_has_ortholog_coords.bed" \
  -v -sorted > "$OUTDIR/human_no_ortholog.humanCoords.bed"