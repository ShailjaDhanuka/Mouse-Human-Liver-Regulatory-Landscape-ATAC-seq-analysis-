#!/bin/bash

#SBATCH --job-name=homer_motif
#SBATCH --partition=RM-shared
#SBATCH --mem=8000M
#SBATCH --cpus-per-task=4
#SBATCH --time=07:00:00

#Usage: bash (or sbatch) enhancer_vs_promoter.sh \
#   <mouse_atac> \
#   <human_atac> \
#   <mouse_peaks_shared_human> \
#   <human_peaks_shared_mouse> \
#   <unique_mouse_human_closed> \
#   <unique_human_mouse_closed> \
#   <mouse_only_no_ortho> \
#   <human_only_no_ortho> \
#   <output_dir>

MOUSE_SHARED_HUMAN_COORDS="$1"
HUMAN_SHARED_MOUSE_COORDS="$2"
UNIQUE_MOUSE_HUMAN_CLOSED="$3"
UNIQUE_HUMAN_MOUSE_CLOSED="$4"
MOUSE_NO_HUMAN="$5"
HUMAN_NO_MOUSE="$6"
OUT_DIR="$7"

# NOTE! - must have hal environment installed to run and source in config (instructions in README)

# fail safe - checks if module exists as command in subshell made by script -- skips if it doesn't to prevent error
command -v module &>/dev/null && module load anaconda3
command -v module &>/dev/null && module load bedtools

# directing conda source to correct file path
source "$(dirname "$0")/../config.sh"
source "$HOMER_CONDA_SOURCE"
conda activate "$HOMER_ENV"

# checking that all input files are correct / non-empty
if [[ $# -ne 7 ]]; then
  echo "Usage: $0 <mouse_shared_human_coords> <human_shared_mouse_coords> <unique_mouse_human_closed> \
  <unique_human_mouse_closed> <mouse_no_human> <human_no_mouse> <output_dir>"
  exit 1
fi

for f in "$MOUSE_SHARED_HUMAN_COORDS" "$HUMAN_SHARED_MOUSE_COORDS" \
"$UNIQUE_MOUSE_HUMAN_CLOSED" "$UNIQUE_HUMAN_MOUSE_CLOSED" "$MOUSE_NO_HUMAN" "$HUMAN_NO_MOUSE"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: input file missing or empty: $f"
    exit 1
  fi
done

mkdir -p $LOCAL/homerRun
mkdir -p $LOCAL/homerRun/preparsed
mkdir -p $OUT_DIR

echo "Copying input files (not genome files) to local node.."
cp $MOUSE_SHARED_HUMAN_COORDS $LOCAL/homerRun/mouse_shared_human_coords
cp $HUMAN_SHARED_MOUSE_COORDS $LOCAL/homerRun/human_shared_mouse_coords
cp $UNIQUE_MOUSE_HUMAN_CLOSED $LOCAL/homerRun/unique_mouse_human_closed
cp $UNIQUE_HUMAN_MOUSE_CLOSED $LOCAL/homerRun/unique_human_mouse_closed
cp $MOUSE_NO_HUMAN $LOCAL/homerRun/mouse_no_human
cp $HUMAN_NO_MOUSE $LOCAL/homerRun/human_no_mouse
echo "Files copied to local node."

echo "Combining open chromatin mapped to closed chromatin orthologs with open regions with no corresponding orthologs.."
cat $LOCAL/homerRun/unique_mouse_human_closed \
    $LOCAL/homerRun/mouse_no_human \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > $LOCAL/homerRun/mouse_only

cat $LOCAL/homerRun/unique_human_mouse_closed \
    $LOCAL/homerRun/human_no_mouse \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > $LOCAL/homerRun/human_only

echo "Combined bed files created."

cd $LOCAL/homerRun

# function to check output of annotated peaks before splitting into promoters vs enhancers
check_annotation() {
    local file=$1
    local label=$2
    if [[ $(wc -l < "$file") -le 1 ]]; then
        echo "ERROR: annotatePeaks.pl produced no results for $label ($file)"
        exit 1
    fi
    echo "Annotation check passed: $label"
}

# function for splitting annotated peak files into promoters vs enhancers for mouse and human
# checks if resulting bed files are empty (0 promoters or enhancers)

split_peaks() {
    local annotated=$1
    local prefix=$2

    awk -F'\t' 'NR>1 && $8 ~ /promoter-TSS/' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_promoter_peaks.bed"

    awk -F'\t' 'NR>1 && ($8 ~ /Intergenic/ || $8 ~ /Intron/)' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_enhancer_peaks.bed"

    if [[ ! -s "${prefix}_promoter_peaks.bed" ]]; then
      echo "WARNING: No promoter peaks found for ${prefix}"
    fi
    if [[ ! -s "${prefix}_enhancer_peaks.bed" ]]; then
      echo "WARNING: No enhancer peaks found for ${prefix}"
    fi

    echo "Created ${prefix}_promoter_peaks.bed and ${prefix}_enhancer_peaks.bed"

}

# annotating mouse and human atac peaks
echo "Annotating peaks..."
annotatePeaks.pl mouse_shared_human_coords hg38 > annotated_shared_human_peaks.txt
annotatePeaks.pl human_shared_mouse_coords mm10 > annotated_shared_mouse_peaks.txt
annotatePeaks.pl mouse_only mm10 > annotated_unique_mouse_peaks.txt
annotatePeaks.pl human_only hg38 > annotated_unique_human_peaks.txt

# checking annotations aren't empty before splitting peaks
check_annotation annotated_shared_mouse_peaks.txt "shared_open_mouse"
check_annotation annotated_shared_human_peaks.txt "shared_open_human"
check_annotation annotated_unique_mouse_peaks.txt "unique_mouse"
check_annotation annotated_unique_human_peaks.txt "unique_human"

echo "Splitting peaks into promoters and enhancers.."
split_peaks annotated_shared_mouse_peaks.txt shared_mouse
split_peaks annotated_shared_human_peaks.txt shared_human
split_peaks annotated_unique_mouse_peaks.txt unique_mouse
split_peaks annotated_unique_human_peaks.txt unique_human

echo "Copying peak BED files to output directory for downstream counting.."
for prefix in shared_mouse shared_human unique_mouse unique_human; do
  for type in promoter enhancer; do
    bed="${prefix}_${type}_peaks.bed"
    [[ -f "$bed" ]] && cp "$bed" "$OUT_DIR/$bed"
  done
done

# COMPARISON 1: shared (foreground) vs unique (background)
# finds what motifs are enriched at conserved open regions vs species-specific ones
echo "Running motif enrichment: shared vs unique (with background).."

for type in promoter enhancer; do

  # first finding enriched motifs in mouse promoters / enhancers
  fg="shared_mouse_${type}_peaks.bed"
  bg="unique_mouse_${type}_peaks.bed"

  # checking that foreground and background sequences arent empty
  if [[ ! -s "$fg" ]]; then
    echo "Skipping $fg (empty)"

  elif [[ ! -s "$bg" ]]; then
    echo "Skipping $fg (background $bg empty)"

  else
    findMotifsGenome.pl "$fg" mm10 \
        $OUT_DIR/shared_vs_unique_mouse_${type}_motifs/ \
        -size 200 -mask \
        -bg "$bg" -useNewBg \
        -preparsedDir $LOCAL/homerRun/preparsed
  fi

  # now finding enriched motifs in human promoters / enhancers
  fg="shared_human_${type}_peaks.bed"
  bg="unique_human_${type}_peaks.bed"

  # checking that foreground and background arent empty
  if [[ ! -s "$fg" ]]; then
    echo "Skipping $fg (empty)"

  elif [[ ! -s "$bg" ]]; then
    echo "Skipping $fg (background $bg empty)"

  else
    findMotifsGenome.pl "$fg" hg38 \
        $OUT_DIR/shared_vs_unique_human_${type}_motifs/ \
        -size 200 -mask \
        -bg "$bg" -useNewBg \
        -preparsedDir $LOCAL/homerRun/preparsed
  fi

done

# COMPARISON 2: unique promoters (fg) vs unique enhancers (bg)
# Finding what motifs distinguish species-specific promoters from species-specific enhancers
echo "Running motif enrichment: unique promoters vs unique enhancers (with background).."

for species in mouse human; do
  [[ $species == "human" ]] && genome=hg38 || genome=mm10

  fg="unique_${species}_promoter_peaks.bed"
  bg="unique_${species}_enhancer_peaks.bed"

  if [[ ! -s "$fg" ]]; then
    echo "Skipping $fg (empty)"
  elif [[ ! -s "$bg" ]]; then
    echo "Skipping $fg (background $bg empty)"
  else
    findMotifsGenome.pl "$fg" $genome \
        $OUT_DIR/unique_${species}_promoter_motifs/ \
        -size 200 -mask \
        -bg "$bg" -useNewBg \
        -preparsedDir $LOCAL/homerRun/preparsed
  fi

  fg="unique_${species}_enhancer_peaks.bed"
  bg="unique_${species}_promoter_peaks.bed"

  if [[ ! -s "$fg" ]]; then
    echo "Skipping $fg (empty)"
  elif [[ ! -s "$bg" ]]; then
    echo "Skipping $fg (background $bg empty)"
  else
    findMotifsGenome.pl "$fg" $genome \
        $OUT_DIR/unique_${species}_enhancer_motifs/ \
        -size 200 -mask \
        -bg "$bg" -useNewBg \
        -preparsedDir $LOCAL/homerRun/preparsed
  fi

done

echo "HOMER done! Results in $OUT_DIR"
