# Mouse-Human Liver Regulatory Landscape — ATAC-seq Analysis
### 03-713: Bioinformatics Data Integration Practicum | Spring 2026

---

## Overview

This project investigates the conservation of transcriptional regulatory activity between human and mouse liver tissue using open chromatin (ATAC-seq) data. We map regulatory elements across species, classify them as enhancers or promoters, identify regulated biological processes, and discover enriched sequence motifs — all tied together in a single automated pipeline.

**Core research questions:**
- To what extent is transcriptional regulatory activity conserved between human and mouse?
- Do enhancers and promoters differ in their degree of cross-species conservation?
- Does the transcriptional regulatory code differ between species and between element types?
- What biological processes are regulated by shared vs. species-specific elements?

---

## Repository Structure

```
.
├── step1_quality_evaluation/        ← Quality assessment of ATAC-seq datasets
├── step2_cross_species_mapping/     ← Liftover & ortholog identification (HALPER)
├── step3_biological_processes/      ← Gene ontology enrichment (rGREAT)
├── step4_enhancers_promoters/       ← Regulatory element classification
├── step5_motif_analysis/            ← Motif discovery with HOMER
├── step6_automated_pipeline/        ← End-to-end automated pipeline
├── scripts/                         ← Shared helper scripts
├── envs/                            ← Conda environment / module files
└── docs/                            ← Documentation, figures, report drafts
```

---

## Pipeline Steps

### Step 1 — Quality Evaluation
Assess the quality of all human and mouse ATAC-seq datasets. Select the highest-quality dataset for downstream analysis based on peak count, FRiP score, and signal enrichment.

**Tools:** bedtools, samtools, R/Python for QC plots

---

### Step 2 — Cross-Species Mapping
Map open chromatin regions between human and mouse genomes using HALPER. Classify regions as shared (ortholog is open in the other species) or species-specific (ortholog is closed).

**Tools:** [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing), [bedtools](https://bedtools.readthedocs.io/en/latest/)

---

### Step 3 — Biological Process Enrichment
Run GO/pathway enrichment on all open chromatin regions, shared regions, and species-specific regions to identify what biological processes are regulated and whether they are conserved.

**Tools:** [rGREAT](https://github.com/jokergoo/rGREAT)

---

### Step 4 — Enhancer and Promoter Classification
Partition open chromatin regions into likely enhancers and promoters. Compare what fraction of each element type is conserved across species.

**Tools:** bedtools, TSS annotations (GENCODE/RefSeq)

---

### Step 5 — Motif Analysis
Discover over-represented sequence motifs in enhancers, promoters, shared regions, and species-specific regions using HOMER.

**Tools:** [HOMER](http://homer.ucsd.edu/homer/)

---

### Step 6 — Automated Pipeline
A single-command pipeline that runs Steps 2–5 sequentially on any Linux cluster with the required tools installed.

```bash
bash step6_automated_pipeline/run_pipeline.sh \
    --human-peaks <path> \
    --mouse-peaks <path> \
    --human-genome <path/to/hg38.fa> \
    --mouse-genome <path/to/mm10.fa> \
    --hal-file <path/to/alignment.hal> \
    --output-dir <output_directory>
```

---

## Tools & References

| Tool | Purpose | Links |
|------|---------|-------|
| HALPER | Cross-species liftover | [GitHub](https://github.com/pfenninglab/halLiftover-postprocessing) · [Paper](https://pubmed.ncbi.nlm.nih.gov/32407523/) |
| bedtools | Genomic interval operations | [Docs](https://bedtools.readthedocs.io/en/latest/) · [Paper](https://pubmed.ncbi.nlm.nih.gov/20110278/) |
| rGREAT | GO enrichment for genomic regions | [GitHub](https://github.com/jokergoo/rGREAT) · [Paper](https://pubmed.ncbi.nlm.nih.gov/36394265/) |
| HOMER | Motif discovery | [Docs](http://homer.ucsd.edu/homer/) · [Paper](https://pubmed.ncbi.nlm.nih.gov/20513432/) |

---

## Contributors

| Name | GitHub |
|------|--------|
| Shailja Dhanuka | [@ShailjaDhanuka](https://github.com/ShailjaDhanuka) |
| | |
| | |

---

## Course

**03-713: Bioinformatics Data Integration Practicum**
Carnegie Mellon University — Spring 2026
