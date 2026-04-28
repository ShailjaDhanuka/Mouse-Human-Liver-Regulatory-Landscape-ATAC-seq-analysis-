# Expected Output Files by Step

All paths below are relative to the `<output_dir>` you pass to the pipeline.
Use these as a reference if you already have results from a given step and want to skip it.

The automated pipeline checks for exactly these files before skipping any step.
If any are missing or empty, the pipeline will exit with an error.

---

## Step 1 — Ortholog Mapping (halLiftover + HALPER)
**Output folder:** `<output_dir>/mapping/`

| File | Description |
|------|-------------|
| `MouseAtacToHumanOrtho.MouseToHuman.HALPER.narrowPeak.gz` | Mouse ATAC peaks lifted over to human (hg38) coordinates |
| `HumanAtacToMouseOrtho.HumanToMouse.HALPER.narrowPeak.gz` | Human ATAC peaks lifted over to mouse (mm10) coordinates |

---

## Step 2 — Shared / Unique Peak Identification (bedtools)
**Output folder:** `<output_dir>/open_chrom/`

| File | Coordinate space | Description |
|------|-----------------|-------------|
| `mouse_shared_in_human.humanCoords.bed` | hg38 | Mouse peaks whose human ortholog is also open |
| `human_shared_in_mouse.mouseCoords.bed` | mm10 | Human peaks whose mouse ortholog is also open |
| `mouse_open_human_closed.mouseCoords.bed` | mm10 | Mouse-specific OCRs (ortholog exists but is closed in human) |
| `human_open_mouse_closed.humanCoords.bed` | hg38 | Human-specific OCRs (ortholog exists but is closed in mouse) |
| `mouse_no_ortholog.mouseCoords.bed` | mm10 | Mouse peaks with no mappable human ortholog |
| `human_no_ortholog.humanCoords.bed` | hg38 | Human peaks with no mappable mouse ortholog |

---

## Step 3 — Gene Ontology Enrichment (rGREAT)
**Output folder:** `<output_dir>/gene_ontology/`

| File | Description                                                     |
|------|-----------------------------------------------------------------|
| `mouse_shared_in_human_GREAT_GO_BP.csv` | GO:BP enrichment for shared mouse peaks in humans (hg38 coords) |
| `human_shared_in_mouse_GREAT_GO_BP.csv` | GO:BP enrichment for shared human peaks in mice (mm10 coords)   |
| `mouse_open_human_closed_GREAT_GO_BP.csv` | GO:BP enrichment for mouse-specific peaks                       |
| `human_open_mouse_closed_GREAT_GO_BP.csv` | GO:BP enrichment for human-specific peaks                       |
| `mouse_no_ortholog_GREAT_GO_BP.csv` | GO:BP enrichment for mouse peaks with no ortholog               |
| `human_no_ortholog_GREAT_GO_BP.csv` | GO:BP enrichment for human peaks with no ortholog               |
| `liver_mouse_GO_BP_allpeaks.csv` | GO:BP enrichment for all mouse ATAC peaks (unfiltered)          |
| `liver_human_GO_BP_allpeaks.csv` | GO:BP enrichment for all human ATAC peaks (unfiltered)          |

---

## Steps 4 & 5 — Enhancer/Promoter Classification + Motif Analysis (HOMER)
**Output folder:** `<output_dir>/homer/`

Each subdirectory below contains a full HOMER `findMotifsGenome` result.
The pipeline checks for `knownResults.txt` and `homerMotifs.all.motifs` inside each directory to verify the run completed.

| Directory                                 | Peak set | Element type |
|-------------------------------------------|----------|-------------|
| `shared_vs_unique_human_enhancer_motifs/` | Shared human peaks | Enhancers |
| `shared_vs_unique_mouse_enhancer_motifs/` | Shared mouse peaks | Enhancers |
| `unique_human_enhancer_motifs/`           | Human-specific peaks | Enhancers |
| `unique_mouse_enhancer_motifs/`           | Mouse-specific peaks | Enhancers |
| `shared_vs_unique_human_promoter_motifs/` | Shared human peaks | Promoters |
| `shared_vs_unique_mouse_promoter_motifs/` | Shared mouse peaks | Promoters |
| `unique_human_promoter_motifs/`           | Human-specific peaks | Promoters |
| `unique_mouse_promoter_motifs/`           | Mouse-specific peaks | Promoters |

---