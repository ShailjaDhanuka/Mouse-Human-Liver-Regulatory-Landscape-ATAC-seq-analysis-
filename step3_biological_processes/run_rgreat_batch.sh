#!/bin/bash
#SBATCH --job-name=rgreat_6files
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH -A bio230007p
#SBATCH --mem=16G

# NOTE! - must have conda environment with rgreat installed and sourced in config (instructions in README)

# fail safe - checks if module exists as command in subshell made by script -- skips if it doesn't to prevent error
command -v module &>/dev/null && module load anaconda3

# directing conda source to correct file path
unset CONDA_PKGS>DIRS
source "$(dirname "$0")/../config.sh"
source "$RGREAT_CONDA_SOURCE"
conda activate "$RGREAT_ENV"

export RGREAT_LIBPATH

# USAGE
# sbatch run_rgreat_batch.sh \
#   <input directory with all 6 bed files from open chromatin identification> \
#   <output directory>

INDIR="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

export INDIR
export OUTDIR

Rscript - <<'EOF'
.libPaths(Sys.getenv("RGREAT_LIBPATH"))

library(rGREAT)
library(GenomicRanges)

INDIR  <- Sys.getenv("INDIR")
OUTDIR <- Sys.getenv("OUTDIR")

# ── helper ────────────────────────────────────────────────────────────────
read_bed_to_gr <- function(path) {
  cat("  Reading:", path, "\n")
  df <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  GRanges(
    seqnames = df[[1]],
    ranges   = IRanges(start = df[[2]] + 1, end = df[[3]]),
    name     = df[[4]]
  )
}

run_great_and_save <- function(bed_path, genome, label, outdir) {
  cat("\n── ", label, " ──────────────────────────\n")
  cat("  Genome:", genome, "\n")

  gr <- read_bed_to_gr(bed_path)
  cat("  Regions loaded:", length(gr), "\n")

  res <- great(gr,
               gene_sets  = "GO:BP",
               tss_source = genome,
               background = NULL,
               cores      = 4)

  tb       <- getEnrichmentTable(res)
  out_path <- file.path(outdir, paste0(label, "_GREAT_GO_BP.csv"))
  write.csv(tb, out_path, row.names = FALSE)
  cat("  Saved:", out_path, "\n")
  cat("  GO terms returned:", nrow(tb), "\n")

  # preview top 5
  cat("  Top 5 terms:\n")
  print(head(tb[, c("id", "description", "p_adjust")], 5))
}

# ── manifest ──────────────────────────────────────────────────────────────
manifest <- list(
  list(
    path   = file.path(INDIR, "mouse_shared_in_human.humanCoords.bed"),
    genome = "hg38",
    label  = "mouse_shared_in_human"
  ),
  list(
    path   = file.path(INDIR, "mouse_open_human_closed.mouseCoords.bed"),
    genome = "mm10",
    label  = "mouse_open_human_closed"
  ),
  list(
    path   = file.path(INDIR, "human_no_ortholog.humanCoords.bed"),
    genome = "hg38",
    label  = "human_no_ortholog"
  ),
  list(
    path   = file.path(INDIR, "human_shared_in_mouse.mouseCoords.bed"),
    genome = "mm10",
    label  = "human_shared_in_mouse"
  ),
  list(
    path   = file.path(INDIR, "human_open_mouse_closed.humanCoords.bed"),
    genome = "hg38",
    label  = "human_open_mouse_closed"
  ),
  list(
    path   = file.path(INDIR, "mouse_no_ortholog.mouseCoords.bed"),
    genome = "mm10",
    label  = "mouse_no_ortholog"
  )
)

# ── batch run ─────────────────────────────────────────────────────────────
cat("Starting GREAT batch run:", length(manifest), "files\n")
cat("Output directory:", OUTDIR, "\n")

for (entry in manifest) {
  tryCatch({
    run_great_and_save(entry$path, entry$genome, entry$label, OUTDIR)
  }, error = function(e) {
    cat("  ERROR on", entry$label, ":", e$message, "\n")
  })
}

# ── final summary ─────────────────────────────────────────────────────────
cat("\n── Final output check ───────────────────────────────\n")
for (entry in manifest) {
  f      <- file.path(OUTDIR, paste0(entry$label, "_GREAT_GO_BP.csv"))
  status <- ifelse(file.exists(f), "OK     ", "MISSING")
  cat(" [", status, "]", basename(f), "\n")
}
EOF