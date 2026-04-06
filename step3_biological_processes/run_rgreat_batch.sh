#!/bin/bash
#SBATCH --job-name=rgreat_6files
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH -A bio230007p
#SBATCH --mem=16G                         
#SBATCH --output=/ocean/projects/bio230007p/sdhanuka/rgreat_6files.log
#SBATCH --error=/ocean/projects/bio230007p/sdhanuka/rgreat_6files.err

unset CONDA_PKGS_DIRS
source /jet/home/sdhanuka/miniconda3/etc/profile.d/conda.sh
conda activate /ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env

# change this to wherever your 6 bedtools output files are
INDIR="/ocean/projects/bio230007p/wanyuef/project/open_chrom_result"
OUTDIR="/ocean/projects/bio230007p/sdhanuka/gene_ontology/batch_GO"
mkdir -p "$OUTDIR"

Rscript - <<'EOF'
.libPaths("/ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env/lib/R/library")

library(rGREAT)
library(GenomicRanges)

INDIR  <- Sys.getenv("INDIR",  "/ocean/projects/bio230007p/wanyuef/project/open_chrom_result")
OUTDIR <- Sys.getenv("OUTDIR", "/ocean/projects/bio230007p/sdhanuka/gene_ontology/batch_GO")

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
    path   = file.path(INDIR, "mouse_open_human_closed.humanCoords.bed"),
    genome = "hg38",
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
    path   = file.path(INDIR, "human_open_mouse_closed.mouseCoords.bed"),
    genome = "mm10",
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