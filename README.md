# ATACer: Mouse-Human Liver Regulatory Landscape — ATAC-seq Analysis
### 03-713: Bioinformatics Data Integration Practicum | Spring 2026

---

## Overview

For our project, we created a pipeline called ATACer that investigates the conservation of transcriptional regulatory activity between human and mouse liver tissue using open chromatin (ATAC-seq) data. We map regulatory elements across species, classify them as enhancers or promoters, identify their biological function, and discover enriched sequence motifs — all tied together in ATACer.

**Core research questions:**
- To what extent is transcriptional regulatory activity conserved between human and mouse liver tissue?
- Do enhancers and promoters differ in their degree of cross-species conservation?
- Does transcriptional regulatory code differ between species and/or between element types?
- What biological processes are regulated by shared vs. species-specific elements?


<img width="1440" height="1312" alt="image" src="https://github.com/user-attachments/assets/9902f431-e00b-4a8c-a879-fd40726eba9e" />

---

## Repository Structure

```
├── sampleOutputs/                   ← Example files of summarized pipeline outputs
├── step1_quality_evaluation/        ← Quality assessment of ATAC-seq datasets
├── step2_cross_species_mapping/     ← Liftover & ortholog identification (halLiftover/HALPER)
├── step3_biological_processes/      ← Gene ontology enrichment (rGREAT)
├── step4-5_enh_vs_prom_and_motifs   ← Regulatory element classification + TF motif analysis (HOMER)
├── step6_automated_pipeline/        ← End-to-end automated pipeline
├── README.md                        ← Pipeline Details
├── config.sh                        ← Configuration file 
├── expected_fileNames.md            ← All expected file names for each step

```

---

## Pipeline Design

### Step 1 — Quality Evaluation
Quality assessment of ATAC-seq data sets were considered when creating this pipeline. Users should ensure they have high-quality ATAC-seq datasets before using the pipeline. The following ENCODE ATAC-seq QC pipeline is recommended for preprocessing raw ATAC reads. ATACer does not perform this step.

**Tools:** [ENCODE-DCC atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)

---

### Step 2 — Cross-Species Mapping
Map open chromatin regions between human and mouse genomes using halLiftover and HALPER. Classify regions as shared (ortholog is open in the other species) or species-specific (ortholog is closed or no ortholog found).

**Tools:** [halLiftover](https://github.com/ComparativeGenomicsToolkit/hal), [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing), [bedtools](https://bedtools.readthedocs.io/en/latest/)

---

### Step 3 — Biological Process Enrichment
Run GO biological pathway enrichment on all open chromatin regions from mouse and human liver, shared regions, and species-specific regions to identify what biological processes are being regulated and whether they are conserved.

**Tools:** [rGREAT](https://github.com/jokergoo/rGREAT)

---

### Step 4 — Enhancer and Promoter Classification
Partition open chromatin regions into likely enhancers and promoters. Compare what fraction of each element type is present in species and also conserved across species.

**Tools:** [HOMER](http://homer.ucsd.edu/homer/)

---

### Step 5 — Motif Analysis
Discover over-represented sequence motifs in species-specific and shared enhancers and promoters regions using HOMER.

**Tools:** [HOMER](http://homer.ucsd.edu/homer/)

---

### Step 6 — Automated Pipeline
A single-command pipeline that runs Steps 2–5 sequentially and provides a summary of results on any Linux cluster with the required tools installed.

```bash
bash step6_automated_pipeline/automatedPipeline.sh \
    <mouse_atac_peaks> \
    <human_atac_peaks> \
    <cactusAlignment> \
    <output_directory> \
    <output_dir> \
    [--start-step <n>]
    [--end-step <n>]
    [--skip <n,n..>]
```
---

## Pipeline Instructions & Usage
### Step 1: Clone the ATACer repo to your working directory.
To use this pipeline, you will have to clone the repository into any desired directory.
```bash
git clone https://github.com/ShailjaDhanuka/Mouse-Human-Liver-Regulatory-Landscape-ATAC-seq-analysis-.git 
```

### Step 2: Dependencies and Installation Instructions
Ensure all necessary dependencies and environments are installed/created. Then update the config.sh file with file paths to your conda source and specific environments. This includes:

- halLiftover and HALPER - [Installation Instructions](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)
  - Create a conda environment and follow the linked instructions to install HAL in the environment
- rGreat - [Installation Instructions](https://github.com/jokergoo/rGREAT/blob/master/README.md) 
  - Create another conda environment and follow the linked instructions to install rGreat in the environment
- HOMER - [Installation Instructions](http://homer.ucsd.edu/homer/introduction/install.html)
  - Create third conda environment and follow the linked instructions to install HOMER in the environment
  - Further installation methods -- can use bioconda to install homer as well, which may be an easier method if you have the bioconda channel
- bedtools -- [Installation Instructions](https://bedtools.readthedocs.io/en/latest/content/installation.html)

### Step 3: Data Preparation
You need to have mouse and human ATAC-seq narrowPeak.gz files, a Cactus HAL alignment file and your desired output folder prepared. 

If you have completed any previous step and would like to skip them in the pipeline, please follow the structure below to make sure output files are placed correctly. Refer to [expected_outputs.md](expected_outputs.md) for specific file information and naming instructions.
- halLiftover/HALPER output (ortholog mapping)-> `<output_dir>/mapping/`
- BED files (open chrom identification) -> `<output_dir>/open_chrom/`
- GO csv files (rGREAT result) -> `<output_dir>/gene_ontology/`
- HOMER BED and motif files (HOMER result) -> `<output_dir>/homer/`


### Step 4: Run the scripts.
ATACer is able to perform 5 possible analyses. You can designate the exact steps you want performed.
- 1 = ortholog mapping (halLiftover/HALPER)
- 2 = shared/unique peaks (bedtools)
- 3 = gene ontology (rGREAT)
- 4 & 5 = enhancer/promoter classification (HOMER annotatePeaks) & motif analysis (HOMER findMotifsGenome)
- 6 = summary and plots (Python)

To run the full pipeline on slurm, follow the instructions here:
```bash
cd /path/to/this/repo/step6_automated_pipelin
sbatch automatedPipeline.sh <mouse_peaks> <human_peaks> <cactus_file> <output_dir>
```
The default for `start-step` and `end-step` is `1` and `6`, respectively. Here is an example of running the pipeline from step 2 to step 6 while skipping steps 3 & 4. You can replace the parameters with whatever you want in your own workflow.
```bash
sbatch automatedPipeline.sh <mouse_peaks> <human_peaks> <cactus_file> <output_dir> --start-step 2 --end-step 6 --skip 3,4
```
The detailed usage instructions are here. Note - there should be no spaces between listed step numbers when using `--skip`:
```bash
sbatch (or bash) automatedPipeline.sh \
    <mouse_peaks> <human_peaks> <cactus_file> <output_dir> \
    [--start-step N] [--end-step N] [--skip N,N,...]
```
The required parameters are:
```bash
    <mouse_peaks> -> path to mouse ATAC-seq
    <human_peaks> -> path to human ATAC-seq
    <cactus_file> -> path to Cactus file
    <output_dir>  -> path to output directory
```
The optional parameters are:
```bash
    --start-step N    -> Start from step N (default: 1)
    --end-step N      -> Stop after step N (default: 6)
    --skip N,N,...    -> Skip specific steps e.g. --skip 3,4
```
## Summary Output
The summarized output will be in your designated output directory under the summary subdirectory. It has the following structure:
```
├── go_enrichment.png                ← top 5 GO enrichment bar plots
├── peak_counts.png                  ← unique/shared peak bar plots
├── promoter_vs_enhancer.png         ← promoter/enhancer peak bar plots
├── pipeline_summary.txt             ← written summary of results
├── stage_log.txt                    ← stage log from pipeline run 
```
You can also go to the [sample output folder](./sampleOutputs/) for reference of what each file contains.

## Tools & References

| Tool          | Purpose                                                      | Links                                                                                                                                              |
|---------------|--------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| halLiftover   | Cross-species genomic mapping                                | [halLiftover GitHub](https://github.com/ComparativeGenomicsToolkit/hal)  |
| HALPER | halLiftover processing for contiguous orthologs construction | [HALPER Github](https://github.com/pfenninglab/halLiftover-postprocessing)
| bedtools      | Genomic interval operations                                  | [Docs](https://bedtools.readthedocs.io/en/latest/) · [Paper](https://pubmed.ncbi.nlm.nih.gov/20110278/)                                            |
| rGREAT        | GO enrichment for genomic regions                            | [GitHub](https://github.com/jokergoo/rGREAT) · [Paper](https://pubmed.ncbi.nlm.nih.gov/36394265/)                                                  |
| HOMER         | Peak annotation and motif discovery                          | [Docs](http://homer.ucsd.edu/homer/) · [Paper](https://pubmed.ncbi.nlm.nih.gov/20513432/)                                                          |

---

## Limitations
There are a few assumptions and limitations to keep in mind when using ATACer.
- The GO analysis does not take peak signal strength into account.
- The full pipeline will take several hours to run on slurm.
- The pipeline assumes that input ATAC-seq data has been properly preprocessed and is of high quality.

## Contributors

| Name               | GitHub                                              |
|--------------------|-----------------------------------------------------|
| Shailja Dhanuka    | [@ShailjaDhanuka](https://github.com/ShailjaDhanuka) |
| Sophia Turecki     | [@sophiat1101](https://github.com/sophiat1101)      |
| Shreya Balamurugan | [@sbalamur02](https://github.com/sbalamur02)        |
| Wanyue Feng        | [@aquatique-plus](https://github.com/aquatique-plus) 
---
## Citation
Balamurugan, S., Dhanuka, S., Feng, W., & Turecki, S. (2026). ATACer: Mouse-Human Liver Regulatory Landscape - ATAC-seq Analysis. 
[03-713 Final Course Project, Carnegie Mellon University] GitHub. https://github.com/ShailjaDhanuka/Mouse-Human-Liver-Regulatory-Landscape-ATAC-seq-analysis-
## Course Project

**03-713: Bioinformatics Data Integration Practicum**
| Carnegie Mellon University — Spring 2026
