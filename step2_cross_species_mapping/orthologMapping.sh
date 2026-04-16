#!/bin/bash

# halper script requires 2000M with 1 cpu -- for two mappings 4000M

#SBATCH -p RM-shared
#SBATCH --mem=4000M
#SBATCH -t 09:00:00
#SBATCH --cpus-per-task=2

# usage:
# sbatch orthologMappping.sh /
#     <MouseAtac.narrowPeak> \
#     <HumanAtac.narrowPeak> \
#     <12mammalianCactus> \
#     <outputDirectory> \

MOUSE_ATAC="$1"
HUMAN_ATAC="$2"
CACTUS_ALIGN="$3"
OUT_DIR="$4"

# NOTE! - must have hal environment installed and source in config (instructions in READ.ME)

## loading in anaconda3 module and hal envs

# fail safe - checks if module exists as command in subshell made by scrip -- skips if it doesnt to prevent error
command -v module &>/dev/null && module load anaconda3
source "$(dirname "$0")/../config.sh"
source "$HAL_CONDA_SOURCE"
export HALPER_DIR
conda activate "$HAL_ENV"

HALPER_SCRIPT="$HALPER_DIR/halper_map_peak_orthologs.sh"

# making local node working directory for job run
mkdir -p $LOCAL/mappingRun

echo "Copying narrowPeak files to local node..."
cp $MOUSE_ATAC $LOCAL/mappingRun/mouse.NarrowPeak.gz
cp $HUMAN_ATAC $LOCAL/mappingRun/human.NarrowPeak.gz
echo "Done copying files to local node."

# in case job fails, sending any results back to project folder
trap 'echo "Copying results back..."; \
      cp $LOCAL/mappingRun/MouseAtacToHumanOrtho* $OUT_DIR 2>/dev/null; \
      cp $LOCAL/mappingRun/HumanAtacToMouseOrtho* $OUT_DIR 2>/dev/null; \
      echo "Done."' EXIT

# running halLiftover/HALPER script in parallel for mouse->human and human->mouse
bash $HALPER_SCRIPT \
-b $LOCAL/mappingRun/mouse.NarrowPeak.gz \
-o $LOCAL/mappingRun \
-s Mouse \
-t Human \
-n MouseAtacToHumanOrtho \
-c $CACTUS_ALIGN &
PID1=$!

bash $HALPER_SCRIPT \
-b $LOCAL/mappingRun/human.NarrowPeak.gz \
-o $LOCAL/mappingRun \
-s Human \
-t Mouse \
-n HumanAtacToMouseOrtho \
-c $CACTUS_ALIGN &
PID2=$!

wait $PID1 $PID2

echo "Copying intermediate files and outputs back to project folder..."
cp $LOCAL/mappingRun/MouseAtacToHumanOrtho* $OUT_DIR
cp $LOCAL/mappingRun/HumanAtacToMouseOrtho* $OUT_DIR
echo "Done"
