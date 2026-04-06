#!/bin/bash

# SBATCH --job-name=homer_motif
# SBATCH --partition=RM-shared
# SBATCH --mem=8GB
# SBATCH --cpus-per-task=4
# SBATCH --time=06:00:00

#Usage: enhancer_vs_promoter.sh \
#   <mouse_atac> \
#   <human_atac> \
#   <mouse_peaks_shared_human> \
#   <human_peaks_shared_mouse> \
#   <unique_mouse_human_closed> \
#   <unique_human_mouse_closed> \
#   <mouse_only_no_ortho> \
#   <human_only_no_ortho> \
#   <human_genome> \
#   <mouse_genome> \
#   <output_dir>

MOUSE_ATAC="$1"
HUMAN_ATAC="$2"
MOUSE_SHARED_HUMAN_COORDS="$3"
HUMAN_SHARED_MOUSE_COORDS="$4"
UNIQUE_MOUSE_HUMAN_CLOSED="$5"
UNIQUE_HUMAN_MOUSE_CLOSED="$6"
MOUSE_NO_HUMAN="$7"
HUMAN_NO_MOUSE="$8"
HUMAN_GENOME="${9}"
MOUSE_GENOME="${10}"
OUT_DIR="${11}"

module load homer/4.11.0
module load bedtools

# checking that all input files specified
if [[ $# -ne 11 ]]; then
  echo "Usage: $0 <mouse_atac> <human_atac> <mouse_shared_human_coords> <human_shared_mouse_coords> <unique_mouse_human_closed> \
  <unique_human_mouse_closed> <mouse_no_human> <human_no_mouse> <human_genome> <mouse_genome> <output_dir>"
  exit 1
fi

# checking for any empty input files
for f in "$MOUSE_ATAC" "$HUMAN_ATAC" "$MOUSE_SHARED_HUMAN_COORDS" "$HUMAN_SHARED_MOUSE_COORDS" \
"$UNIQUE_MOUSE_HUMAN_CLOSED" "$UNIQUE_HUMAN_MOUSE_CLOSED" "$MOUSE_NO_HUMAN" "$HUMAN_NO_MOUSE"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: input file missing or empty: $f"
    exit 1
  fi
done

# making local directories to hold input files and outputs during job
mkdir -p $LOCAL/homerRun
mkdir -p $LOCAL/homerRun/preparsed
mkdir -p $OUT_DIR

echo "Copying input files (not genome files) to local node.."
cp $MOUSE_ATAC $LOCAL/homerRun/mouse_atac
cp $HUMAN_ATAC $LOCAL/homerRun/human_atac
cp $MOUSE_SHARED_HUMAN_COORDS $LOCAL/homerRun/mouse_shared_human_coords
cp $HUMAN_SHARED_MOUSE_COORDS $LOCAL/homerRun/human_shared_mouse_coords
cp $UNIQUE_MOUSE_HUMAN_CLOSED $LOCAL/homerRun/unique_mouse_human_closed
cp $UNIQUE_HUMAN_MOUSE_CLOSED $LOCAL/homerRun/unique_human_mouse_closed
cp $MOUSE_NO_HUMAN $LOCAL/homerRun/mouse_no_human
cp $HUMAN_NO_MOUSE $LOCAL/homerRun/human_no_mouse
echo "Files copied to local node."

# combining the two open chromatin files that are unique to only one species
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
    local threshold=${3:-1000}

    awk -F'\t' -v thresh="$threshold" \
        'NR>1 && ($10>=-thresh && $10<=thresh)' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_promoter_peaks.bed"

    if [[ ! -s "${prefix}_promoter_peaks.bed" ]]; then
      echo "WARNING: No promoter peaks found for ${prefix}"
    fi

    awk -F'\t' -v thresh="$threshold" \
        'NR>1 && ($10<-thresh || $10>thresh)' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_enhancer_peaks.bed"

    if [[ ! -s "${prefix}_enhancer_peaks.bed" ]]; then
      echo "WARNING: No enhancer peaks found for ${prefix}"
    fi

    echo "Created ${prefix}_promoter_peaks.bed and ${prefix}_enhancer_peaks.bed with threshold ${threshold}bp"

}

# annotating mouse and human atac peaks
echo "Annotating peaks..."
annotatePeaks.pl mouse_atac $MOUSE_GENOME > annotated_mouse_peaks.txt
annotatePeaks.pl human_atac $HUMAN_GENOME > annotated_human_peaks.txt
annotatePeaks.pl mouse_shared_human_coords $HUMAN_GENOME > annotated_shared_human_peaks.txt
annotatePeaks.pl human_shared_mouse_coords $MOUSE_GENOME > annotated_shared_mouse_peaks.txt
annotatePeaks.pl mouse_only $MOUSE_GENOME > annotated_unique_mouse_peaks.txt
annotatePeaks.pl human_only $HUMAN_GENOME > annotated_unique_human_peaks.txt

# checking annotations aren't empty before splitting peaks
check_annotation annotated_mouse_peaks.txt "mouse_atac"
check_annotation annotated_human_peaks.txt "human_atac"
check_annotation annotated_shared_mouse_peaks.txt "shared_open_mouse"
check_annotation annotated_shared_human_peaks.txt "shared_open_human"
check_annotation annotated_unique_mouse_peaks.txt "unique_mouse"
check_annotation annotated_unique_human_peaks.txt "unique_human"

echo "Splitting peaks into promoters and enhancers.."
split_peaks annotated_mouse_peaks.txt mouse 1000
split_peaks annotated_human_peaks.txt human 1000
split_peaks annotated_shared_mouse_peaks.txt shared_mouse 1000
split_peaks annotated_shared_human_peaks.txt shared_human 1000
split_peaks annotated_unique_mouse_peaks.txt unique_mouse 1000
split_peaks annotated_unique_human_peaks.txt unique_human 1000

echo "Running motif enrichment.."

# looping over each prefix of files to do motif analysis
for prefix in mouse human shared_mouse shared_human unique_mouse unique_human; do

  # fail safe -- if bed file is empty, skipping motif analysis
  for type in promoter enhancer; do
    bed="${prefix}_${type}_peaks.bed"
    if [[ ! -s "$bed" ]]; then
      echo "Skipping motif analysis for $bed (empty)"
      continue
    fi

    if [[ $prefix == "human" || $prefix == "shared_human" || $prefix == "unique_human" ]]; then
        genome=$HUMAN_GENOME
    else
        genome=$MOUSE_GENOME
    fi

    # calling motif enrichment function in homer centered 200bp around peak summits
    # masking repetitive/low complexity sequences
    findMotifsGenome.pl "$bed" $genome \
        $OUT_DIR/${prefix}_${type}_motifs/ \
        -size 200 -mask \
        -preparsedDir $LOCAL/homerRun/preparsed

  done
done

echo "Done! Results in $OUT_DIR"
