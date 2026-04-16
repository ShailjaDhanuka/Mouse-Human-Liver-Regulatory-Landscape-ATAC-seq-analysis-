#!/bin/bash
#SBATCH --job-name=rgreat_analysis
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH -A bio230007p

module load anaconda3

unset CONDA_PKGS_DIRS
source "$(dirname "$0")/../config.sh"
source "$RGREAT_CONDA_SOURCE"
conda activate "$RGREAT_ENV"

export RGREAT_LIBPATH

MOUSE_ATAC="$1"
HUMAN_ATAC="$2"
OUTPUT_DIR="$3"

export MOUSE_ATAC
export HUMAN_ATAC
export OUTPUT_DIR

Rscript - <<'EOF'
.libPaths(Sys.getenv("RGREAT_LIBPATH"))

library(rGREAT)
library(GenomicRanges)

#  ---- HUMAN ----
cat("Running GREAT on human peaks...\n")
human_peaks <- read.table(
  gzfile(Sys.getenv("HUMAN_ATAC")),
  header = FALSE,
  col.names = c("chr","start","end","name","score","strand",
                 "fold_enrichment","pval","qval","summit")
 )

human_gr <- GRanges(
  seqnames = human_peaks$chr,
  ranges   = IRanges(start = human_peaks$start, end = human_peaks$end)
 )

human_res <- great(human_gr, "GO:BP", "hg38")
human_go  <- getEnrichmentTable(human_res)
write.csv(human_go,
  file.path(Sys.getenv("OUTPUT_DIR"), "liver_human_GO_BP_allpeaks.csv"), row.names = FALSE)
cat("Human done!\n")

# ---- MOUSE ----
cat("Running GREAT on mouse peaks...\n")
mouse_peaks <- read.table(
  gzfile(Sys.getenv("MOUSE_ATAC")),
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
  file.path(Sys.getenv("OUTPUT_DIR"), "liver_mouse_GO_BP_allpeaks.csv"),
  row.names = FALSE)
cat("Mouse done!\n")

# ---- PREVIEW ----
cat("\nTop human GO terms:\n")
print(head(human_go[, c("id","description","p_adjust")]))
cat("\nTop mouse GO terms:\n")
print(head(mouse_go[, c("id","description","p_adjust")]))
EOF
