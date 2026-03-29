#!/bin/bash
#SBATCH --job-name=rgreat_analysis
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH -A bio230007p
#SBATCH --output=/ocean/projects/bio230007p/sdhanuka/rgreat.log
#SBATCH --error=/ocean/projects/bio230007p/sdhanuka/rgreat.err



#this is where  the script to execute rgreat is (i ran it from my /ocean/projects/bio230007p/sdhanuka/ as a batch script)
unset CONDA_PKGS_DIRS
source /jet/home/sdhanuka/miniconda3/etc/profile.d/conda.sh
conda activate /ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env

Rscript - <<'EOF'
.libPaths("/ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env/lib/R/library")

library(rGREAT)

# # ---- HUMAN ----
# cat("Running GREAT on human peaks...\n")
# human_peaks <- read.table(
#   gzfile("/ocean/projects/bio230007p/ikaplow/HumanAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"),
#   header = FALSE,
#   col.names = c("chr","start","end","name","score","strand",
#                 "fold_enrichment","pval","qval","summit")
# )

# human_gr <- GRanges(
#   seqnames = human_peaks$chr,
#   ranges   = IRanges(start = human_peaks$start, end = human_peaks$end)
# )

# human_res <- great(human_gr, "GO:BP", "hg38")
# human_go  <- getEnrichmentTable(human_res)
# write.csv(human_go,
#   "/ocean/projects/bio230007p/sdhanuka/gene_ontology/liver_human_GO_BP_allpeaks.csv",
#   row.names = FALSE)
# cat("Human done!\n")

# ---- MOUSE ----
cat("Running GREAT on mouse peaks...\n")
mouse_peaks <- read.table(
  gzfile("/ocean/projects/bio230007p/ikaplow/MouseAtac/Liver/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz"),
  header = FALSE,
  col.names = c("chr","start","end","name","score","strand",
                "fold_enrichment","pval","qval","summit")
)

mouse_gr <- GRanges(
  seqnames = mouse_peaks$chr,
  ranges   = IRanges(start = mouse_peaks$start, end = mouse_peaks$end)
)

mouse_res <- great(mouse_gr, "GO:BP", "mm10")
mouse_go  <- getEnrichmentTable(mouse_res)
write.csv(mouse_go,
  "/ocean/projects/bio230007p/sdhanuka/gene_ontology/liver_mouse_GO_BP_allpeaks.csv",
  row.names = FALSE)
cat("Mouse done!\n")

# ---- PREVIEW ----
cat("\nTop human GO terms:\n")
print(head(human_go[, c("id","description","p_adjust")]))
cat("\nTop mouse GO terms:\n")
print(head(mouse_go[, c("id","description","p_adjust")]))
EOF
