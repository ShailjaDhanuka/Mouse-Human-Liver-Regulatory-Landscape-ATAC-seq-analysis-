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

# NOTE! -- HALPER_DIR must be set to path of halLiftover-processing repo
# example: export HALPER_DIR=/ocean/projects/bio23000p/USER/halLiftover-postprocessing

if [[ -z "$HALPER_DIR" ]]; then
      echo "Error: HALPER_DIR is not set in script." >&2
      echo "Set it to path of halLiftover-postprocessing repo." >&2
      echo "Ex. export HALPER_DIR=/path/to/halLiftover-postprocessing?" >&2
      exit 1
fi

HALPER_SCRIPT="$HALPER_DIR/halper_map_peak_orthologs.sh"

# NOTE! - must have hal environment installed from halLiftover-postprocessing instructions
# or create env from this repo (ex. conda env create -f halEnv.yml)

# loading in anaconda3 module and hal envs
module load anaconda3
conda activate hal

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
