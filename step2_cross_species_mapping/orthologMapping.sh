#!/bin/bash

# halper script requires 2000M with 1 cpu -- for two mappings 4000M

#SBATCH -p RM-shared
#SBATCH --mem=4000M
#SBATCH -t 09:00:00
#SBATCH --cpus-per-task=2

# loading in anaconda3 module  and hal envs
module load anaconda3
source activate hal

# making local working directory for job run
mkdir -p $LOCAL/mappingRun

echo "Copying narrowPeak files to local node..."
cp /ocean/projects/bio230007p/ikaplow/MouseAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz $LOCAL/mappingRun/mouse.NarrowPeak.gz
cp /ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz $LOCAL/mappingRun/human.NarrowPeak.gz
echo "Done copying files to local node."

# in case job fails, sending any results back to project folder
trap 'echo "Copying results back..."; \
      cp $LOCAL/mappingRun/MouseAtacToHumanOrtho* /ocean/projects/bio230007p/sturecki/project/mappingResults/ 2>/dev/null; \
      cp $LOCAL/mappingRun/HumanAtacToMouseOrtho* /ocean/projects/bio230007p/sturecki/project/mappingResults/ 2>/dev/null; \
      echo "Done."' EXIT

# running halLiftover/HALPER script in parallel for mouse->human and human->mouse
bash /jet/home/sturecki/practicum/projectScripts/halper_map_peak_orthologs.sh \
-b $LOCAL/mappingRun/mouse.NarrowPeak.gz \
-o $LOCAL/mappingRun \
-s Mouse \
-t Human \
-n MouseAtacToHumanOrtho \
-c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal &
PID1=$!

bash /jet/home/sturecki/practicum/projectScripts/halper_map_peak_orthologs.sh \
-b $LOCAL/mappingRun/human.NarrowPeak.gz \
-o $LOCAL/mappingRun \
-s Human \
-t Mouse \
-n HumanAtacToMouseOrtho \
-c /ocean/projects/bio230007p/ikaplow/Alignments/10plusway-master.hal &
PID2=$!

wait $PID1 $PID2

echo "Copying intermediate and outputs back to project folder..."
cp $LOCAL/mappingRun/MouseAtacToHumanOrtho* /ocean/projects/bio230007p/sturecki/project/mappingResults/
cp $LOCAL/mappingRun/HumanAtacToMouseOrtho* /ocean/projects/bio230007p/sturecki/project/mappingResults/
echo "Done"
