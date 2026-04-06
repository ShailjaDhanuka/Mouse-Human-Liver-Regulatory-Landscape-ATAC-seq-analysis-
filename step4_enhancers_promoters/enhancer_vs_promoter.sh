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
#   <output_dir>

MOUSE_ATAC="$1"
HUMAN_ATAC="$2"
MOUSE_SHARED_HUMAN_COORDS="$3"
HUMAN_SHARED_MOUSE_COORDS="$4"
UNIQUE_MOUSE_HUMAN_CLOSED="$5"
UNIQUE_HUMAN_MOUSE_CLOSED="$6"
MOUSE_NO_HUMAN="$7"
HUMAN_NO_MOUSE="$8"
OUT_DIR="$9"

module load homer
module load bedtools

# checking that all input parameters were entered
if [[ $# -ne 9 ]]; then
  echo "Usage: $0 <mouse_atac> <human_atac> <mouse_shared_human_coords> <human_shared_mouse_coords> <unique_mouse_human_closed> \
  <unique_human_mouse_closed> <mouse_no_human> <human_no_mouse> <output_dir>"
  exit 1
fi

# checking that all files contain data 
for f in "$MOUSE_ATAC" "$HUMAN_ATAC" "$MOUSE_SHARED_HUMAN_COORDS" "$HUMAN_SHARED_MOUSE_COORDS" \
"$UNIQUE_MOUSE_HUMAN_CLOSED" "$UNIQUE_HUMAN_MOUSE_CLOSED" "$MOUSE_NO_HUMAN" "$HUMAN_NO_MOUSE"; do
  if [[ ! -s "$f" ]]; then
    echo "ERROR: input file missing or empty: $f"
    exit 1
  fi
done

# making local directory to hold intermediate outputs and making output directory if doesn't exist yet
mkdir -p $LOCAL/homerRun
mkdir -p $OUT_DIR

# copying files to local node for easy access during job run
echo "Copying input files to local node.."
cp $MOUSE_ATAC $LOCAL/homerRun/mouse_atac
cp $HUMAN_ATAC $LOCAL/homerRun/human_atac
cp $MOUSE_SHARED_HUMAN_COORDS $LOCAL/homerRun/mouse_shared_human_coords
cp $HUMAN_SHARED_MOUSE_COORDS $LOCAL/homerRun/human_shared_mouse_coords
cp $UNIQUE_MOUSE_HUMAN_CLOSED $LOCAL/homerRun/unique_mouse_human_closed
cp $UNIQUE_HUMAN_MOUSE_CLOSED $LOCAL/homerRun/unique_human_mouse_closed
cp $MOUSE_NO_HUMAN $LOCAL/homerRun/mouse_no_human
cp $HUMAN_NO_MOUSE $LOCAL/homerRun/human_no_mouse
echo "Files copied to local node."

# combining the files for open chromatin with closed orthologs or open chromatin with no orthologs for both species to get
# one file of unique open chromatin each
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
    # checking if file is not empty after annotation
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

    # splitting tab deliminated columns and checking if 10th col (distance to TSS) is less than 1000bp away from peak
    awk -F'\t' -v thresh="$threshold" \
        'NR>1 && ($10>=-thresh && $10<=thresh)' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_promoter_peaks.bed"

    # checking if promoters were found
    if [[ ! -s "${prefix}_promoter_peaks.bed" ]]; then
      echo "WARNING: No promoter peaks found for ${prefix}"
    fi
    
    # same - splitting tab deliminated columns but checking if 10th col (distance to TSS) is more than 1000bp away from peak
    awk -F'\t' -v thresh="$threshold" \
        'NR>1 && ($10<-thresh || $10>thresh)' "$annotated" | \
        awk -F'\t' 'OFS="\t" {print $2, $3, $4}' > "${prefix}_enhancer_peaks.bed"

    # checking if enhancers were found
    if [[ ! -s "${prefix}_enhancer_peaks.bed" ]]; then
      echo "WARNING: No enhancer peaks found for ${prefix}"
    fi

    echo "Created ${prefix}_promoter_peaks.bed and ${prefix}_enhancer_peaks.bed with threshold ${threshold}bp"

}

# annotating mouse and human atac peaks with HOMER
echo "Annotating peaks..."
annotatePeaks.pl mouse_atac mm39 > annotated_mouse_peaks.txt
annotatePeaks.pl human_atac hg38 > annotated_human_peaks.txt
annotatePeaks.pl mouse_shared_human_coords hg38 > annotated_shared_human_peaks.txt
annotatePeaks.pl human_shared_mouse_coords mm39 > annotated_shared_mouse_peaks.txt
annotatePeaks.pl mouse_only mm39 > annotated_unique_mouse_peaks.txt
annotatePeaks.pl human_only hg38 > annotated_unique_human_peaks.txt

# checking annotations aren't empty before splitting peaks
check_annotation annotated_mouse_peaks.txt "mouse_atac"
check_annotation annotated_human_peaks.txt "human_atac"
check_annotation annotated_shared_mouse_peaks.txt "shared_open_mouse"
check_annotation annotated_shared_human_peaks.txt "shared_open_human"
check_annotation annotated_unique_mouse_peaks.txt "unique_mouse"
check_annotation annotated_unique_human_peaks.txt "unique_human"

# splitting annotated peaks into promoters or enhancers with HOMER
echo "Splitting peaks into promoters and enhancers.."
split_peaks annotated_mouse_peaks.txt mouse 1000
split_peaks annotated_human_peaks.txt human 1000
split_peaks annotated_shared_mouse_peaks.txt shared_mouse 1000
split_peaks annotated_shared_human_peaks.txt shared_human 1000
split_peaks annotated_unique_mouse_peaks.txt unique_mouse 1000
split_peaks annotated_unique_human_peaks.txt unique_human 1000


# moving onto HOMER motif enrichment analysis for all files made after splitting
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

    # getting right genome code
    if [[ $prefix == "human" || $prefix == "shared_human" || $prefix == "unique_human"]]; then
        genome=hg38
    else
        genome=mm39
    fi

    # calling motif enrichment function in homer centered 200bp around peak summits
    # masking repetitive/low complexity sequences
    findMotifsGenome.pl "$bed" $genome \
        $OUT_DIR/${prefix}_${type}_motifs/ \
        -size 200 -mask

  done
done

echo "Done! Results in $OUT_DIR"