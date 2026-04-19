# ATACer: Mouse-Human Liver Regulatory Landscape — ATAC-seq Analysis
### 03-713: Bioinformatics Data Integration Practicum | Spring 2026

---

## Overview

This project creates a pipeline called ATACer and investigates the conservation of transcriptional regulatory activity between human and mouse liver tissue using open chromatin (ATAC-seq) data. We map regulatory elements across species, classify them as enhancers or promoters, identify their biological function, and discover enriched sequence motifs — all tied together in ATACer.

**Core research questions:**
- To what extent is transcriptional regulatory activity conserved between human and mouse?
- Do enhancers and promoters differ in their degree of cross-species conservation?
- Does the transcriptional regulatory code differ between species and between element types?
- What biological processes are regulated by shared vs. species-specific elements?


<img width="1440" height="1312" alt="image" src="https://github.com/user-attachments/assets/9902f431-e00b-4a8c-a879-fd40726eba9e" />

---

## Repository Structure

```
├── step1_quality_evaluation         ← Quality assessment of ATAC-seq datasets
├── step1_cross_species_mapping/     ← Liftover & ortholog identification (HALPER)
├── step2_biological_processes/      ← Gene ontology enrichment (rGREAT)
├── step3-4_enh_vs_prom_and_motifs   ← Regulatory element classification + TF motif analysis (HOMER)
├── step5_automated_pipeline/        ← End-to-end automated pipeline
├── README.md                        ← Pipeline Details
├── config.sh                        ← Configuration file 

```

---

## Pipeline Design

### Step 1 — Quality Evaluation
Quality assessment outputs for ATAC-seq data sets considered for creation of this pipeline. Select the highest-quality dataset for downstream analysis based on peak count, FRiP score, and signal enrichment.

**Tools:** [ENCODE-DCC atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)

---

### Step 2 — Cross-Species Mapping
Map open chromatin regions between human and mouse genomes using halLiftover and HALPER. Classify regions as shared (ortholog is open in the other species) or species-specific (ortholog is closed).

**Tools:** [halLiftover](https://github.com/ComparativeGenomicsToolkit/hal), [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing), [bedtools](https://bedtools.readthedocs.io/en/latest/)

---

### Step 3 — Biological Process Enrichment
Run GO/pathway enrichment on all open chromatin regions, shared regions, and species-specific regions to identify what biological processes are regulated and whether they are conserved.

**Tools:** [rGREAT](https://github.com/jokergoo/rGREAT)

---

### Step 4 — Enhancer and Promoter Classification
Partition open chromatin regions into likely enhancers and promoters. Compare what fraction of each element type is conserved across species.

**Tools:** [HOMER](http://homer.ucsd.edu/homer/)

---

### Step 5 — Motif Analysis
Discover over-represented sequence motifs in enhancers, promoters, shared regions, and species-specific regions using HOMER.

**Tools:** [HOMER](http://homer.ucsd.edu/homer/)

---

### Step 6 — Automated Pipeline
A single-command pipeline that runs Steps 2–5 sequentially on any Linux cluster with the required tools installed.

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
### Step 1: Clone this repo to your working folder.
To begin using this pipeline, you will have to clone the repository into you desired folder.
```bash
git clone https://github.com/ShailjaDhanuka/Mouse-Human-Liver-Regulatory-Landscape-ATAC-seq-analysis-.git 
```

### Step 2: Dependencies and Installation Instructions
Ensure all necessary dependencies and environments are installed/created.  Then update config.sh file with file paths to your conda source and specific environments. This includes:

- halLiftover and HALPER (install in its own hal environment) - [Installation Instructions](https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/hal_install_instructions.md)
- rGreat (install in its own rGreat environment) - [Installation Instructions](https://github.com/jokergoo/rGREAT/blob/master/README.md)
- HOMER (install in its own HOMER environment) - easiest to use bioconda to install (must have bioconda channel)
  - ``conda install bioconda::homer``
  - Further installation methods - [Installation Instructions](http://homer.ucsd.edu/homer/introduction/install.html)
- bedtools -- [Installation Instructions](https://bedtools.readthedocs.io/en/latest/content/installation.html)

### Step 3: Prepare your data
You need to have mouse&humane ATAC-seq files, cactus files and your desired output folder prepared. If you have done mapping, open chromatin identification or GO enrichment, please follow the instructions here to make sure output files are placed correctly:
- HALPER output (ortholog mapping)-> `<output_dir>/mapping/`
- BED files (open chrom identification) -> `<output_dir>/open_chrom/`
- GO csv files (rGREAT result) -> `<output_dir>/gene_ontology/`


### Step 4: Run the scripts.
This pipeline is able to perform 6 possible downstream analysis. And you can designate the exact steps you want.
- 1 = ortholog mapping (HALPER)
- 2 = shared/unique peaks (bedtools)
- 3 = gene ontology (rGREAT)
- 4 = enhancer/promoter classification (HOMER annotatePeaks)
- 5 = motif analysis (HOMER findMotifsGenome)
- 6 = summary and plots (Python)

To run the pipeline on slurm, follow the intsructions here:
```bash
cd /path/to/this/respo/
cd ./step6_automated_pipelin
sbatch automatedPipeline.sh <mouse_peaks> <human_peaks> <cactus_file> <output_dir> --start-step 1 --end-step 6 --skip 3,4
```
The default set for `start-step` and `end-step` is `1` and `6` separately. Here is an example of running the pipeline from step 1 to step 6 while skipping step 3&4. You can replace the parameters with whatever you want in your own workflow.

The detailed usage instructions are here:
```bash
sbatch (or bash) automatedPipeline.sh \
    <mouse_peaks> <human_peaks> <cactus_file> <output_dir> \
    [--start-step N] [--end-step N] [--skip N,N,...]
```
The required parameters are:
```bash
    <mouse_peaks> -> path to mouse ATAC-seq
    <human_peaks> -> path to human ATAC-seq
    <cactus_file> -> path to cactus file
    <output_dir>  -> path to output directory
```
The optional parameters are:
```bash
    --start-step N    -> Start from step N (default: 1)
    --end-step N      -> Stop after step N (default: 6)
    --skip N,N,...    -> Skip specific steps e.g. --skip 3,4
```
## Sample Output
The output should be in your designated output directory in the following structure.
```
├── go_enrichment.png                ← GO enrichment results
├── peak_counts.png                  ← unique/shared peak counting
├── promoter_vs_enhancer.png         ← promoter/enhancer peak counting
├── pipeline_summary.txt             ← summary of all the results
├── stage_log.txt                    ← log file
```
You can also go to the [sample output folder](./sampleOutputs/) for reference.

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
When using this pipeline, there are a few assumptions and limitations to keep in mind.
- The GO analysis does not take peak signal strength into account.
- The pipeline assumes that the input data has been properly preprocessed and is of high quality.

## Contributors

| Name               | GitHub                                              |
|--------------------|-----------------------------------------------------|
| Shailja Dhanuka    | [@ShailjaDhanuka](https://github.com/ShailjaDhanuka) |
| Sophia Turecki     | [@sophiat1101](https://github.com/sophiat1101)      |
| Shreya Balamurugan | [@sbalamur02](https://github.com/sbalamur02)        |
| Wanyue Feng        | [@aquatique-plus](https://github.com/aquatique-plus) 
---

## Course Project

**03-713: Bioinformatics Data Integration Practicum**
| Carnegie Mellon University — Spring 2026
