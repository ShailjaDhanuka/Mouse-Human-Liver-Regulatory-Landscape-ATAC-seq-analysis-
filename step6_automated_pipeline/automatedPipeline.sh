#!/bin/bash

# stop the pipeline if any command fails, if an undefined variable is used,
# or if a command inside a pipe fails
set -euo pipefail
# make sure the trap catches errors inside functions too
set -E

# this pipeline compares open chromatin regions between human and mouse liver
# it runs 6 steps in order: ortholog mapping, peak comparison, gene ontology,
# enhancer/promoter classification, motif analysis, and a final summary

# slurm job settings - tells bridges2 how much memory and time we need
#SBATCH --job-name=liver_pipeline
#SBATCH -p RM-shared
#SBATCH -t 15:00:00
#SBATCH --mem=12000M
#SBATCH --cpus-per-task=6
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# ============================================================
# HOW TO RUN
# ============================================================
# run all steps:
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder>
#
# run only steps 1 and 2:
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --end-step 2
#
# skip step 3 (rGREAT):
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --skip 3
#
# start from step 2 (ortholog mapping already done):
#   NOTE: HALPER narrowPeak files must be in <output_folder>/mapping/
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --start-step 2
#
# start from step 3 (bedtools already done):
#   NOTE: BED files must be in <output_folder>/open_chrom/
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --start-step 3
#
# start from step 4 (rGREAT already done):
#   NOTE: BED files must be in <output_folder>/open_chrom/
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --start-step 4
#
# start from step 6 (HOMER already done):
#   NOTE: GO CSV files must be in <output_folder>/gene_ontology/
#   AND HOMER result directories must be in <output_folder>/homer/
# sbatch automatedPipeline.sh \
#     <mouse_peaks> <human_peaks> <cactus_file> <output_folder> \
#     --start-step 6

# ============================================================
# CHECK INPUTS
# ============================================================

# make sure the user gave us at least 4 arguments
# if not, print a help message and exit
if [ "$#" -lt 4 ]; then
    echo "ERROR: Not enough arguments"
    echo "Usage: $0 <mouseAtac> <humanAtac> <cactusFile> <outputFolder>"
    echo ""
    echo "Optional flags:"
    echo "  --start-step N    Start from step N (default: 1)"
    echo "  --end-step N      Stop after step N (default: 6)"
    echo "  --skip N,N,...    Skip specific steps e.g. --skip 3,4"
    echo ""
    echo "Steps:"
    echo "  1 = ortholog mapping (HALPER)"
    echo "  2 = shared/unique peaks (bedtools)"
    echo "  3 = gene ontology (rGREAT)"
    echo "  4 = enhancer/promoter classification (HOMER annotatePeaks)"
    echo "  5 = motif analysis (HOMER findMotifsGenome)"
    echo "  6 = summary and plots (Python)"
    exit 1
fi

# save each input argument as a named variable so they are easy to use later
MOUSE_ATAC="$1"    # mouse atac-seq peaks file
HUMAN_ATAC="$2"    # human atac-seq peaks file
CACTUS_FILE="$3"   # cactus alignment hal file for cross-species mapping
OUTPUT_DIR="$4"    # folder where all results will be saved
shift 4            # move past the 4 required arguments so we can read optional flags

# ============================================================
# STEP CONTROL
# ============================================================

# default values - run all steps unless the user says otherwise
START_STEP=1
END_STEP=6
SKIP_STEPS=""

# read any optional flags the user provided
# loops through remaining arguments and reads --start-step, --end-step, --skip
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --start-step)
            START_STEP="$2"   # start from this step number
            shift 2
            ;;
        --end-step)
            END_STEP="$2"     # stop after this step number
            shift 2
            ;;
        --skip)
            SKIP_STEPS="$2"   # comma separated list of steps to skip
            shift 2
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            echo "Valid options: --start-step, --end-step, --skip"
            exit 1
            ;;
    esac
done

# ============================================================
# SETUP
# ============================================================

# create all output folders upfront so each step has somewhere to save results
# -p means dont error if the folder already exists
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/mapping"        # step 1 outputs go here
mkdir -p "$OUTPUT_DIR/open_chrom"     # step 2 outputs go here
mkdir -p "$OUTPUT_DIR/gene_ontology"  # step 3 outputs go here
mkdir -p "$OUTPUT_DIR/homer"          # steps 4 and 5 outputs go here
mkdir -p "$OUTPUT_DIR/summary"        # step 6 outputs go here

# path to the stage log file - tracks which steps passed, failed or were skipped
STAGE_LOG="$OUTPUT_DIR/summary/stage_log.txt"

# write pipeline run info to the stage log
echo "Pipeline started: $(date)" > "$STAGE_LOG"
echo "Mouse ATAC:  $MOUSE_ATAC"  >> "$STAGE_LOG"
echo "Human ATAC:  $HUMAN_ATAC"  >> "$STAGE_LOG"
echo "Cactus file: $CACTUS_FILE" >> "$STAGE_LOG"
echo "Output dir:  $OUTPUT_DIR"  >> "$STAGE_LOG"
echo ""                          >> "$STAGE_LOG"

# load the tools we need on bridges2
module load anaconda3
module load bedtools

# ============================================================
# TRAP
# ============================================================

# this variable tracks which step we are currently running
# it gets updated before each step so if the pipeline crashes
# we know exactly where it failed
CURRENT_STEP="starting up"

# this function runs automatically if the pipeline crashes
# it prints where the failure happened, lists any files saved
# so far, and tells the user how to resume from where it stopped
cleanup_on_error() {
    local exit_code=$?
    echo ""
    echo "pipeline failed at: $CURRENT_STEP (exit code $exit_code)"
    echo ""
    echo "files saved so far:"

    # check each output folder and print how many files are in it
    for dir in mapping open_chrom gene_ontology homer summary; do
        if [ -d "$OUTPUT_DIR/$dir" ]; then
            local count
            count=$(find "$OUTPUT_DIR/$dir" -type f | wc -l)
            echo "  $OUTPUT_DIR/$dir ($count files)"
        fi
    done

    # write the failure info to the stage log
    echo "PIPELINE FAILED at: $CURRENT_STEP (exit code: $exit_code)" >> "$STAGE_LOG"
    echo "failed at: $(date)" >> "$STAGE_LOG"

    echo ""
    echo "to resume, run:"
    echo "  sbatch automatedPipeline.sh $MOUSE_ATAC $HUMAN_ATAC $CACTUS_FILE $OUTPUT_DIR --start-step $START_STEP"
    echo "see README for required files per step"
}

# register the trap - cleanup_on_error will run automatically if any command fails
trap cleanup_on_error ERR

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# returns true if a step should run, false if it should be skipped
# checks start step, end step, and skip list
should_run() {
    local step=$1

    # skip if before the start step
    if [ "$step" -lt "$START_STEP" ]; then return 1; fi

    # skip if after the end step
    if [ "$step" -gt "$END_STEP" ]; then return 1; fi

    # skip if in the skip list
    if [[ -n "$SKIP_STEPS" ]]; then
        local IFS=','
        read -ra SKIPPED <<< "$SKIP_STEPS"
        for s in "${SKIPPED[@]}"; do
            if [ "$step" -eq "$s" ]; then return 1; fi
        done
    fi

    # if none of the above conditions matched, the step should run
    return 0
}

# writes a line to the stage log recording whether a step passed,
# failed, or was skipped, along with the current timestamp
log_stage() {
    local step="$1"
    local status="$2"    # PASSED, FAILED, or SKIPPED
    local message="$3"
    echo "[$status] Step $step: $message ($(date))" >> "$STAGE_LOG"
}

# checks that a list of output files exist and are not empty
# if any file is missing, prints an error and exits the pipeline
# also checks gzip files are not corrupted
check_outputs() {
    local description="$1"   # plain english description of what we are checking
    shift
    local files=("$@")       # list of files we expect to exist

    for f in "${files[@]}"; do
        # check the file exists and is not empty (-s means exists and size > 0)
        if [[ ! -s "$f" ]]; then
            echo "ERROR: $description"
            echo "  missing or empty: $f"
            echo "  if skipping a step, make sure its output files"
            echo "  are in the correct folder. see README for details."
            exit 1
        fi

        # if the file is gzip compressed, check it is not corrupted
        if [[ "$f" == *.gz ]]; then
            gzip -t "$f" || {
                echo "ERROR: corrupted gzip file: $f"
                exit 1
            }
        fi
    done
    echo "check passed: $description"
}

# checks that a BED file is valid - has the right number of columns,
# enough lines, and valid numeric coordinates
assert_valid_bed() {
    local file="$1"
    local description="$2"
    local min_lines="${3:-1}"   # minimum number of lines expected (default 1)

    # check file exists and is not empty
    if [[ ! -s "$file" ]]; then
        echo "ERROR: $description - file missing or empty: $file"
        exit 1
    fi

    # check file has at least 3 columns (chromosome, start, end)
    local cols
    cols=$(awk '{print NF; exit}' "$file")
    if [[ "$cols" -lt 3 ]]; then
        echo "ERROR: $description - not a valid BED file"
        echo "  expected at least 3 columns, found $cols in: $file"
        exit 1
    fi

    # check file has at least the minimum number of lines (peaks)
    local lines
    lines=$(wc -l < "$file")
    if [[ "$lines" -lt "$min_lines" ]]; then
        echo "ERROR: $description - too few peaks"
        echo "  expected at least $min_lines lines, found $lines in: $file"
        exit 1
    fi

    # check that start and end coordinates are valid numbers
    if ! awk '{if($2!~/^[0-9]+$/ || $3!~/^[0-9]+$/){exit 1}}' "$file"; then
        echo "ERROR: $description - coordinates are not numbers in: $file"
        exit 1
    fi

    # check that start coordinate is always less than end coordinate
    if awk '{if($2>=$3){exit 1}}' "$file"; then
        echo "ERROR: $description - some peaks have start >= end in: $file"
        exit 1
    fi

    echo "check passed: $file looks like a valid BED file ($lines peaks)"
}

# ============================================================
# PRINT RUN PLAN
# ============================================================

# print a summary of what will run and what will be skipped
echo "-------------------------------------------"
echo "liver regulatory pipeline"
echo "-------------------------------------------"
echo "mouse atac: $MOUSE_ATAC"
echo "human atac: $HUMAN_ATAC"
echo "output:     $OUTPUT_DIR"
echo ""

# names for each step used in the run plan printout and stage log
STEP_NAMES=("" 
    "ortholog mapping (HALPER)"
    "shared/unique peaks (bedtools)"
    "gene ontology (rGREAT)"
    "enhancer/promoter classification (HOMER annotatePeaks)"
    "motif analysis (HOMER findMotifsGenome)"
    "summary and plots (Python)"
)

# loop through each step and print whether it will run or be skipped
for i in 1 2 3 4 5 6; do
    if should_run $i; then
        echo "  step $i: ${STEP_NAMES[$i]} -- will run"
        log_stage $i "PLANNED" "${STEP_NAMES[$i]}"
    else
        echo "  step $i: ${STEP_NAMES[$i]} -- skipping"
        log_stage $i "SKIPPED" "${STEP_NAMES[$i]}"
    fi
done

# print a reminder if any steps are being skipped
if [ "$START_STEP" -gt 1 ] || [ -n "$SKIP_STEPS" ]; then
    echo ""
    echo "note: some steps are being skipped."
    echo "      make sure required input files are already"
    echo "      in the correct output folder."
    echo "      see README for details."
fi

echo "-------------------------------------------"

# ============================================================
# STEP 1 - ortholog mapping (HALPER)
# ============================================================

if should_run 1; then

    # update current step so the trap knows where we are if something fails
    CURRENT_STEP="step 1 - ortholog mapping (HALPER)"
    echo ""
    echo "--- step 1: ortholog mapping (HALPER) ---"

    # check all three input files exist before trying to run
    for f in "$MOUSE_ATAC" "$HUMAN_ATAC" "$CACTUS_FILE"; do
        if [[ ! -s "$f" ]]; then
            echo "ERROR: input file missing or empty: $f"
            exit 1
        fi
    done

    # run the HALPER mapping script in both directions
    # mouse to human and human to mouse
    bash orthologMapping.sh \
        "$MOUSE_ATAC" \
        "$HUMAN_ATAC" \
        "$CACTUS_FILE" \
        "$OUTPUT_DIR/mapping"

    # verify both HALPER output files were produced and are not empty
    check_outputs \
        "HALPER narrowPeak files (step 1 - ortholog mapping output) in $OUTPUT_DIR/mapping/" \
        "$OUTPUT_DIR/mapping/MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz" \
        "$OUTPUT_DIR/mapping/HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz"

    # count how many peaks were successfully mapped and print to screen
    MOUSE_MAPPED=$(zcat "$OUTPUT_DIR/mapping/MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz" | wc -l)
    HUMAN_MAPPED=$(zcat "$OUTPUT_DIR/mapping/HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz" | wc -l)
    echo "  mouse peaks mapped to human: $MOUSE_MAPPED"
    echo "  human peaks mapped to mouse: $HUMAN_MAPPED"

    log_stage 1 "PASSED" "ortholog mapping done. mouse mapped: $MOUSE_MAPPED human mapped: $HUMAN_MAPPED"
    echo "step 1 done"
else
    echo "skipping step 1 - ortholog mapping (HALPER)"

    # if step 1 is skipped but later steps need its outputs,
    # check that the HALPER files are already present in the mapping folder
    if should_run 2 || should_run 3 || should_run 4 || should_run 5; then
        check_outputs \
            "HALPER narrowPeak files (step 1 - ortholog mapping output) must be in $OUTPUT_DIR/mapping/ when skipping step 1" \
            "$OUTPUT_DIR/mapping/MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz" \
            "$OUTPUT_DIR/mapping/HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz"
    fi
fi

# ============================================================
# STEP 2 - shared/unique peaks (bedtools)
# ============================================================

if should_run 2; then
    CURRENT_STEP="step 2 - shared/unique peaks (bedtools)"
    echo ""
    echo "--- step 2: shared/unique peaks (bedtools) ---"

    # check that HALPER files from step 1 are present before running
    # these are needed as input whether step 1 ran now or previously
    check_outputs \
        "HALPER narrowPeak files (step 1 - ortholog mapping output) required in $OUTPUT_DIR/mapping/" \
        "$OUTPUT_DIR/mapping/MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz" \
        "$OUTPUT_DIR/mapping/HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz"

    # run bedtools to compare peaks and split into shared and unique categories
    bash open_chromatin_finding_script.sh \
        "$OUTPUT_DIR/mapping/MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz" \
        "$OUTPUT_DIR/mapping/HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz" \
        "$MOUSE_ATAC" \
        "$HUMAN_ATAC" \
        "$OUTPUT_DIR/open_chrom"

    # verify all 6 output BED files were produced
    check_outputs \
        "BED files (step 2 - shared/unique peaks output) in $OUTPUT_DIR/open_chrom/" \
        "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_no_ortholog.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_no_ortholog.humanCoords.bed"

    # validate the format of the 4 main BED files
    # (skipping no_ortholog files since they may have dots in column 4)
    for bed in \
        "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed"; do
        assert_valid_bed "$bed" "BED file (step 2 - shared/unique peaks output)" 1
    done

    # count peaks in each category and print a summary
    MOUSE_SHARED=$(wc -l < "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed")
    HUMAN_SHARED=$(wc -l < "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed")
    MOUSE_UNIQUE=$(wc -l < "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed")
    HUMAN_UNIQUE=$(wc -l < "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed")
    MOUSE_NO_ORTHO=$(wc -l < "$OUTPUT_DIR/open_chrom/mouse_no_ortholog.mouseCoords.bed")
    HUMAN_NO_ORTHO=$(wc -l < "$OUTPUT_DIR/open_chrom/human_no_ortholog.humanCoords.bed")

    echo "  mouse shared in human:    $MOUSE_SHARED peaks"
    echo "  human shared in mouse:    $HUMAN_SHARED peaks"
    echo "  mouse open human closed:  $MOUSE_UNIQUE peaks"
    echo "  human open mouse closed:  $HUMAN_UNIQUE peaks"
    echo "  mouse no ortholog:        $MOUSE_NO_ORTHO peaks"
    echo "  human no ortholog:        $HUMAN_NO_ORTHO peaks"

    log_stage 2 "PASSED" "shared/unique peaks done. mouse shared: $MOUSE_SHARED human shared: $HUMAN_SHARED mouse unique: $MOUSE_UNIQUE human unique: $HUMAN_UNIQUE"
    echo "step 2 done"
else
    echo "skipping step 2 - shared/unique peaks (bedtools)"

    # if step 2 is skipped but later steps need its outputs,
    # check that the BED files are already present in the open_chrom folder
    if should_run 3 || should_run 4 || should_run 5; then
        check_outputs \
            "BED files (step 2 - shared/unique peaks output) must be in $OUTPUT_DIR/open_chrom/ when skipping step 2" \
            "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
            "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed" \
            "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed" \
            "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed" \
            "$OUTPUT_DIR/open_chrom/mouse_no_ortholog.mouseCoords.bed" \
            "$OUTPUT_DIR/open_chrom/human_no_ortholog.humanCoords.bed"
    fi
fi

# ============================================================
# STEP 3 - gene ontology (rGREAT)
# ============================================================

if should_run 3; then
    CURRENT_STEP="step 3 - gene ontology (rGREAT)"
    echo ""
    echo "--- step 3: gene ontology (rGREAT) ---"

    # check that BED files from step 2 are present before running
    check_outputs \
        "BED files (step 2 - shared/unique peaks output) required in $OUTPUT_DIR/open_chrom/" \
        "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed"

    # run rGREAT on all 6 BED files to find enriched biological processes
    bash run_rgreat_batch.sh \
        "$OUTPUT_DIR/open_chrom" \
        "$OUTPUT_DIR/gene_ontology"

    # verify all 6 GO result CSV files were produced
    check_outputs \
        "GO result CSV files (step 3 - gene ontology output) in $OUTPUT_DIR/gene_ontology/" \
        "$OUTPUT_DIR/gene_ontology/mouse_shared_in_human_GREAT_GO_BP.csv" \
        "$OUTPUT_DIR/gene_ontology/human_shared_in_mouse_GREAT_GO_BP.csv" \
        "$OUTPUT_DIR/gene_ontology/mouse_open_human_closed_GREAT_GO_BP.csv" \
        "$OUTPUT_DIR/gene_ontology/human_open_mouse_closed_GREAT_GO_BP.csv" \
        "$OUTPUT_DIR/gene_ontology/mouse_no_ortholog_GREAT_GO_BP.csv" \
        "$OUTPUT_DIR/gene_ontology/human_no_ortholog_GREAT_GO_BP.csv"

    log_stage 3 "PASSED" "gene ontology done. 6 GO result files produced"
    echo "step 3 done"
else
    echo "skipping step 3 - gene ontology (rGREAT)"

    # if step 3 is skipped but step 6 needs its outputs,
    # check that the GO CSV files are already present
    if should_run 6; then
        check_outputs \
            "GO result CSV files (step 3 - gene ontology output) must be in $OUTPUT_DIR/gene_ontology/ when skipping step 3" \
            "$OUTPUT_DIR/gene_ontology/mouse_shared_in_human_GREAT_GO_BP.csv" \
            "$OUTPUT_DIR/gene_ontology/human_shared_in_mouse_GREAT_GO_BP.csv" \
            "$OUTPUT_DIR/gene_ontology/mouse_open_human_closed_GREAT_GO_BP.csv" \
            "$OUTPUT_DIR/gene_ontology/human_open_mouse_closed_GREAT_GO_BP.csv" \
            "$OUTPUT_DIR/gene_ontology/mouse_no_ortholog_GREAT_GO_BP.csv" \
            "$OUTPUT_DIR/gene_ontology/human_no_ortholog_GREAT_GO_BP.csv"
    fi
fi

# ============================================================
# STEPS 4 AND 5 - enhancer/promoter and motif analysis (HOMER)
# ============================================================

if should_run 4 || should_run 5; then
    CURRENT_STEP="steps 4/5 - enhancer/promoter and motif analysis (HOMER)"
    echo ""
    echo "--- steps 4/5: enhancer/promoter classification and motif analysis (HOMER) ---"

    # check that BED files from step 2 are present before running
    check_outputs \
        "BED files (step 2 - shared/unique peaks output) required in $OUTPUT_DIR/open_chrom/" \
        "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_no_ortholog.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_no_ortholog.humanCoords.bed"

    # run HOMER to classify peaks as enhancers or promoters
    # and find enriched transcription factor motifs in each set
    bash enhancer_vs_promoter.sh \
        "$OUTPUT_DIR/open_chrom/mouse_shared_in_human.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_shared_in_mouse.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_open_human_closed.humanCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_open_mouse_closed.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/mouse_no_ortholog.mouseCoords.bed" \
        "$OUTPUT_DIR/open_chrom/human_no_ortholog.humanCoords.bed" \
        "$OUTPUT_DIR/homer"

    # verify all 8 HOMER motif result directories were produced
    # (4 peak categories x 2 element types = 8 directories)
    check_outputs \
        "HOMER motif result files (steps 4/5 - enhancer/promoter and motif analysis output) in $OUTPUT_DIR/homer/" \
        "$OUTPUT_DIR/homer/shared_human_enhancer_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/shared_mouse_enhancer_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/unique_human_enhancer_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/unique_mouse_enhancer_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/shared_human_promoter_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/shared_mouse_promoter_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/unique_human_promoter_motifs/knownResults.txt" \
        "$OUTPUT_DIR/homer/unique_mouse_promoter_motifs/knownResults.txt"

    log_stage 4 "PASSED" "enhancer/promoter classification done"
    log_stage 5 "PASSED" "motif analysis done"
    echo "steps 4 and 5 done"
else
    echo "skipping steps 4 and 5 - enhancer/promoter classification and motif analysis (HOMER)"

    # if steps 4/5 are skipped but step 6 needs their outputs,
    # check that the HOMER result files are already present
    if should_run 6; then
        check_outputs \
            "HOMER motif result files (steps 4/5 - enhancer/promoter and motif analysis output) must be in $OUTPUT_DIR/homer/ when skipping steps 4/5" \
            "$OUTPUT_DIR/homer/shared_human_enhancer_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/shared_mouse_enhancer_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/unique_human_enhancer_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/unique_mouse_enhancer_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/shared_human_promoter_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/shared_mouse_promoter_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/unique_human_promoter_motifs/knownResults.txt" \
            "$OUTPUT_DIR/homer/unique_mouse_promoter_motifs/knownResults.txt"
    fi
fi

# ============================================================
# STEP 6 - summary and plots (Python)
# ============================================================

if should_run 6; then
    CURRENT_STEP="step 6 - summary and plots (Python)"
    echo ""
    echo "--- step 6: summary and plots (Python) ---"

    # run the python summary script which reads all results and
    # generates a text summary and bar charts
    python summarize_results.py \
        --openChromDir "$OUTPUT_DIR/open_chrom" \
        --goDir        "$OUTPUT_DIR/gene_ontology" \
        --homerDir     "$OUTPUT_DIR/homer" \
        --outputDir    "$OUTPUT_DIR/summary"

    # verify the summary text file and peak count plot were produced
    check_outputs \
        "summary files (step 6 - summary and plots output) in $OUTPUT_DIR/summary/" \
        "$OUTPUT_DIR/summary/pipeline_summary.txt" \
        "$OUTPUT_DIR/summary/peak_counts.png"

    log_stage 6 "PASSED" "summary and plots done"
    echo "step 6 done"
else
    echo "skipping step 6 - summary and plots (Python)"
fi

# ============================================================
# DONE
# ============================================================

# write the finish time to the stage log
echo "" >> "$STAGE_LOG"
echo "pipeline finished: $(date)" >> "$STAGE_LOG"

# print a final summary to the screen
echo ""
echo "-------------------------------------------"
echo "pipeline done"
echo "-------------------------------------------"
echo ""
echo "stage log:"
cat "$STAGE_LOG"
echo ""
echo "results saved to:"
echo "  $OUTPUT_DIR/mapping/        HALPER narrowPeak files"
echo "  $OUTPUT_DIR/open_chrom/     BED files (shared and unique peaks)"
echo "  $OUTPUT_DIR/gene_ontology/  GO result CSV files"
echo "  $OUTPUT_DIR/homer/          HOMER motif result directories"
echo "  $OUTPUT_DIR/summary/        summary text, plots and stage log"
echo "-------------------------------------------"