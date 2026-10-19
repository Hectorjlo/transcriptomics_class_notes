# Transcriptomics Analysis & Computational Workflows

**Author:** Héctor Jesús López-Ordaz  
**Institution:** Licenciatura en Ciencias Genómicas (LCG), Universidad Nacional Autónoma de México (UNAM)

---

## Overview

This repository contains end-to-end bioinformatic workflows, benchmarking analyses, and reproducible scripts developed for bulk and single-cell RNA sequencing (RNA-seq) data analysis. 

The core project (`DE_analysis/`) investigates age-related transcriptional remodeling in murine bone tissue (**3-month-old vs. 24-month-old male *Mus musculus***) using data from the *Tabula Muris Senis* bulk RNA-seq dataset ([**GSE132040**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132040), *"Ageing hallmarks exhibit organ-specific temporal signatures"*). Beyond the biological characterization of skeletal aging and *inflammaging*, this repository implements a systematic **12-way methodological benchmark** evaluating how sequencing architecture (*Single-End* vs. *Paired-End*), read mapping paradigms (*HISAT2*, *STAR*, and *Salmon*), and statistical frameworks (*DESeq2* vs. *edgeR*) influence differential gene expression discovery.

Additionally, the repository includes workflow automation modules in **Nextflow (DSL2)**, a custom alignment-log parser written in **Rust**, multifactorial differential expression designs, and **single-cell RNA-seq (scRNA-seq)** clustering and cell-type annotation pipelines using **Seurat** and **SingleR**.

---

## Repository Structure

```text
transcriptomics_class_notes/
├── DE_analysis/                  # Core Project: Murine Bone Aging (GSE132040, 3m vs. 24m)
│   ├── data/                     # Metadata (SRA accessions), GENCODE vM36 maps, exclusion BEDs, DAVID outputs
│   ├── src/                      # Modular pipeline scripts (Bash, Python, Rust, R)
│   │   ├── extractSRR.py         # Automated NCBI GEO metadata filtering for target SRR accessions
│   │   ├── auto_fastp.sh         # Parallelized adapter trimming and quality filtering (fastp)
│   │   ├── auto_hisat2*.sh       # Splice-aware genome alignment in SE and PE modes (HISAT2)
│   │   ├── auto_star*.sh         # Ultrafast splice-aware genome alignment in SE and PE modes (STAR)
│   │   ├── auto_BAMic.sh         # Parallel SAM-to-BAM conversion, MAPQ >= 10 filtering, and sorting (samtools)
│   │   ├── auto_pileup_gen.sh    # BPM-normalized BigWig coverage track generation with blacklist filtering (deepTools)
│   │   ├── auto_ft_counts.sh     # Gene-level read quantification (Subread featureCounts)
│   │   ├── auto_salmon.sh        # Transcriptome pseudoalignment in SE and PE modes (Salmon)
│   │   ├── autrge.R              # GFF3 parsing and transcript-to-gene summarization (tximport)
│   │   ├── combine.R             # Aggregation of sample-level featureCounts outputs into global count matrices
│   │   ├── DEA_DESeq2.R          # CLI-driven Differential Expression Analysis pipeline using DESeq2
│   │   ├── DEA_EdgeR.R           # CLI-driven Differential Expression Analysis pipeline using edgeR (GLM LRT)
│   │   ├── run_all.sh            # Master runner executing all 12 Aligner × Mode × Statistical Model combinations
│   │   ├── GO_enrichment.R       # Gene Ontology (BP/CC) bubble plot visualization from DAVID outputs
│   │   ├── GSEA.R                # Ranked Gene Set Enrichment Analysis via clusterProfiler & enrichplot
│   │   ├── plot_stats.R          # Alignment benchmark visualization
│   │   └── stats/                # Rust CLI utility for parsing and normalizing STAR & HISAT2 summary logs
│   ├── results/                  # Output matrices, DEG tables, BigWig tracks, and publication-grade figures
│   ├── Tarea1_LopezOrdaz_HectorJesus.md   # Report 1: QC, Trimming, and Genome Alignment Benchmark
│   ├── Tarea2_LopezOrdaz_HectorJesus.qmd  # Report 2: BAM Processing, Pileups, featureCounts & Pseudoalignment
│   ├── Tarea4_LopezOrdaz_HectorJesus.qmd  # Report 3: Comparative Differential Expression (DESeq2 vs. edgeR)
│   └── Tarea5_LopezOrdaz_HectorJesus.qmd  # Report 4: Functional Enrichment (DAVID, GSEA) & PPI Networks (STRING)
│
├── nextflow_examples/            # Reproducible Workflow Management (Nextflow DSL2)
│   ├── src/fastpc.nf             # Parallelized raw QC (FastQC), read trimming (fastp), and post-trim QC subworkflows
│   └── results/                  # Nextflow execution reports, MultiQC summaries, and status heatmaps
│
└── DGE_ex/                       # Specialized Differential Expression & Single-Cell RNA-seq Case Studies
    ├── c8/                       # Bulk RNA-seq DEA workflow (4m vs. 18m across sexes)
    ├── c9/                       # Multifactorial experimental design with interaction terms (~ genotype * treatment)
    └── c12/                      # Single-cell RNA-seq pipeline (Seurat, SingleR automated annotation, MAST/Wilcoxon/DESeq2)
```

---

## Core Project: Transcriptomic Profiling of Murine Bone Aging (`DE_analysis/`)

### 1. Experimental Design & Preprocessing
Seven male bone RNA-seq libraries (`SRR9126859`, `SRR9126934`, `SRR9127063`, `SRR9127293`, `SRR9127308`, `SRR9127357`, `SRR9127396`; Illumina NovaSeq 6000, 100 bp) were programmaticallly selected from **GSE132040** (`src/extractSRR.py`), comparing young adult (**3 months**, $n=3$) and aged (**24 months**, $n=4$) cohorts.

* **Quality Control & Trimming:** Raw reads were profiled with `FastQC` and `MultiQC`, revealing 5' technical bias (first 10–15 bp), *Nextera* adapter contamination, and *Clontech SMARTer II A* oligonucleotides. Parallelized processing with `fastp` (`src/auto_fastp.sh`) using sliding-window quality trimming ($Q \ge 20$), poly-G tail removal, and adapter detection restored per-sequence quality metrics across all samples.

| Raw Read Quality (`MultiQC`) | Post-Trimming Read Quality (`MultiQC`) |
| :---: | :---: |
| ![Raw FastQC Status](DE_analysis/results/plots/raw_fastqc-status-check-heatmap.png) | ![Clean FastQC Status](DE_analysis/results/plots/clean_fastqc-status-check-heatmap.png) |

---

### 2. Alignment, Pseudoalignment & Computational Benchmarking
To evaluate the trade-offs between traditional splice-aware genome aligners and lightweight transcriptome pseudoaligners, reads were mapped against the **GRCm39 (mm39 / GENCODE vM36)** reference in both **Single-End (SE)** and **Paired-End (PE)** modes:
1. **Genome Alignment:** `STAR` and `HISAT2` (allowing up to 7 multi-mappings), followed by parallel `samtools` filtering ($\text{MAPQ} \ge 10$), coordinate sorting, indexing (`src/auto_BAMic.sh`), BPM-normalized `BigWig` track generation excluding `mm39` problematic regions (`deepTools bamCoverage`, `src/auto_pileup_gen.sh`), and gene-level counting (`Subread featureCounts`, `src/auto_ft_counts.sh`).
2. **Transcriptome Pseudoalignment:** `Salmon` (`src/auto_salmon.sh`) coupled with custom GFF3 parsing and transcript-to-gene summarization via `tximport` (`src/autrge.R`).
3. **Log Normalization in Rust:** A custom Rust CLI tool (`src/stats/`) parses heterogeneous log files from `STAR` and `HISAT2` into a unified record format for downstream visualization in R (`src/plot_stats.R`).

| Step / Modality | Avg. CPU Threads | Wall-Clock Time | Key Observation |
| :--- | :---: | :---: | :--- |
| **`fastp` (Trimming)** | ~18.3 | ~3 min | High parallel efficiency via GNU `parallel` |
| **`STAR` (Single-End)** | ~6.9 | ~35 min | Consistently higher unique mapping rate (+5–10% vs. `HISAT2`) |
| **`STAR` (Paired-End)** | ~6.0 | ~1 h 49 min | Optimal balance of speed and spatial concordance |
| **`HISAT2` (Single-End)** | ~3.9 | ~1 h 38 min | Lower memory footprint; stricter unaligned rates |
| **`HISAT2` (Paired-End)** | ~2.5 | ~15 h 17 min | High computational cost under strict PE concordance |

![Alignment Statistics Comparison](DE_analysis/results/plots/stats_points.png)

---

### 3. 12-Way Differential Expression Benchmark (`DESeq2` vs. `edgeR`)
Command-line R workflows (`src/DEA_DESeq2.R` and `src/DEA_EdgeR.R`) were engineered using `optparse` to run identically structured analyses across all **12 combinations** ($\{\text{HISAT2}, \text{STAR}, \text{Salmon}\} \times \{\text{SE}, \text{PE}\} \times \{\text{DESeq2}, \text{edgeR}\}$), applying $|\log_2\text{FC}| > 0.5$ and $\text{FDR} < 0.05$:

| Aligner / Quantifier | Statistical Method | Read Architecture | Upregulated DEGs | Downregulated DEGs | Total DEGs |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **HISAT2** | `DESeq2` | Single-End | 48 | 6 | 54 |
| **HISAT2** | `DESeq2` | Paired-End | 23 | 7 | 30 |
| **HISAT2** | `edgeR` | Single-End | 56 | 8 | 64 |
| **HISAT2** | `edgeR` | Paired-End | 20 | 9 | 29 |
| **STAR** | `DESeq2` | Single-End | 47 | 6 | 53 |
| **STAR** | `DESeq2` | Paired-End | **36** | **8** | **44** |
| **STAR** | `edgeR` | Single-End | 63 | 8 | 71 |
| **STAR** | `edgeR` | Paired-End | 48 | 8 | 56 |
| **Salmon** (`tximport`) | `DESeq2` | Single-End | 47 | 10 | 57 |
| **Salmon** (`tximport`) | `DESeq2` | Paired-End | 38 | 10 | 48 |
| **Salmon** (`tximport`) | `edgeR` | Single-End | 36 | 7 | 43 |
| **Salmon** (`tximport`) | `edgeR` | Paired-End | 29 | 6 | 35 |

#### Methodological Takeaways
* **Paired-End vs. Single-End Precision:** Single-end mode consistently inflated DEG counts due to unconstrained read placement. Paired-end mode enforced mate orientation and insert-size concordance, reducing false-positive calls in repetitive or homologous loci.
* **Algorithmic Specialization:** `edgeR` (GLM Likelihood Ratio Test) exhibited higher sensitivity on genome-aligned count matrices (`STAR`/`HISAT2`), whereas `DESeq2` captured more DEGs when integrated with `Salmon`/`tximport` offsets, effectively accounting for transcript-level abundance uncertainty.
* **Core Biological Reproducibility:** Despite pipeline variations, key extracellular matrix markers (*Acan*, *Col10a1*, *Col9a3*, *Chad*), signaling regulators (*Rasd2*), and immunoglobulin transcripts (*Ighg2c*) were reproducibly identified across workflows.

| `STAR` + `DESeq2` (Paired-End) Top 20 DEGs | `STAR` + `DESeq2` (Paired-End) Volcano Plot |
| :---: | :---: |
| ![Top 20 Heatmap](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/plots/heatmap_top20.png) | ![Volcano Plot](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/plots/vulcano_plot.png) |

---

### 4. Functional Enrichment, GSEA & Protein Interactome
Using the high-confidence `STAR` + `DESeq2` (Paired-End) gene set, functional characterization was performed via **DAVID** (Gene Ontology Biological Process & Cellular Component), **clusterProfiler** (Ranked GSEA), and **STRING** (Protein-Protein Interaction networks):

1. **Extracellular Matrix Collapse vs. Humoral *Inflammaging* (GO BP & CC):**
   * **Downregulated programs** are strongly enriched for *cartilage development*, *skeletal system development*, *collagen trimer*, and *extracellular matrix* structural components.
   * **Upregulated programs** show massive enrichment ($\text{Fold Enrichment} > 200$) for *circulating immunoglobulin complexes*, *complement activation*, and *humoral immune response*, reflecting age-associated chronic inflammation (*inflammaging*).
2. **Osteoclast Remodeling Dynamics (GSEA):**
   * Threshold-free GSEA revealed simultaneous positive enrichment of **osteoclast proliferation** ($\text{NES} = 1.832, p_{\text{adj}} < 0.01$) and **negative regulation of osteoclast differentiation** ($\text{NES} = 1.600, p_{\text{adj}} < 0.05$), suggesting an expansion of osteoclast progenitors paired with compensatory feedback attempting to restrain bone resorption.
3. **Interactome Topology (STRING PPI):**
   * Downregulated genes form a tightly interconnected structural module (*Col10a1*–*Col9a3*–*Acan*–*Chad*–*Cnmd*), demonstrating coordinated shutdown of the cartilage/bone matrix machinery, whereas upregulated effectors form a dispersed secretome highlighted by the *Cxcl13*–*Ccl8* chemokine axis.

| GO Biological Process (Upregulated) | GO Biological Process (Downregulated) |
| :---: | :---: |
| ![GO Up BP](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/GO_Bubble_Plot_upregulated_BP.png) | ![GO Down BP](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/GO_Bubble_Plot_downregulated_BP.png) |

| Ranked GSEA (Osteoclast Pathways) | STRING PPI Network (Combined DEGs) |
| :---: | :---: |
| ![GSEA Plot](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/GSEA_plot.png) | ![STRING PPI](DE_analysis/results/star/DE_analysis/DESeq2/paired_end/Functional_analysis/plots/STRING_DEG.png) |

---

## Additional Modules

### Nextflow DSL2 Preprocessing Pipeline (`nextflow_examples/`)
[`nextflow_examples/src/fastpc.nf`](nextflow_examples/src/fastpc.nf) encapsulates raw read quality assessment (`FastQC`), adapter/quality trimming (`fastp`), and post-trimming validation into reusable Nextflow DSL2 subworkflows (`raw_qc`, `trim_by_fastp`, `clean_qc`) using `channel.fromFilePairs` for concurrent execution.

### Multifactorial & Single-Cell RNA-seq Workflows (`DGE_ex/`)
* **Complex Experimental Designs (`DGE_ex/c9/`):** Implements generalized linear models with interaction terms (`~ genotype + treatment + genotype:treatment`) in `DESeq2` to disentangle genotype-specific treatment responses (WT vs. KO).
* **Single-Cell RNA-seq Pipeline (`DGE_ex/c12/scRNAseq.R`):** End-to-end 10x Genomics scRNA-seq analysis using `Seurat` (mitochondrial/doublet QC filtering, `vst` variable feature selection, PCA, Louvain clustering, UMAP embedding), automated reference-based cell type annotation via `SingleR` (`celldex` mouse RNA-seq reference), and comparative marker detection between Microglia and Neurons using **Wilcoxon Rank-Sum**, **MAST**, and **DESeq2**.

---

## Technologies & Dependencies

* **Languages:** R, Bash (GNU Parallel), Python 3, Rust (2021 Edition), Nextflow (DSL2), Quarto / LaTeX
* **Bioinformatics CLI Tools:** `sra-tools`, `FastQC`, `MultiQC`, `fastp`, `STAR`, `HISAT2`, `Salmon`, `samtools`, `deepTools` (`bamCoverage`), `Subread` (`featureCounts`)
* **R / Bioconductor Ecosystem:** `DESeq2`, `edgeR`, `tximport`, `ComplexHeatmap`, `clusterProfiler`, `enrichplot`, `org.Mm.eg.db`, `Seurat`, `SingleR`, `celldex`, `presto`, `ggplot2`, `dplyr`, `optparse`

---

## Reproducibility & Usage

To execute the complete 12-way differential expression benchmark once count tables and `tximport` objects are generated:

```bash
cd DE_analysis
chmod +x src/run_all.sh
./src/run_all.sh
```

To compile the Rust alignment-summary parser and generate normalized alignment metrics:

```bash
cd DE_analysis/src/stats
cargo build --release
./target/release/stats ../../results/ > ../../results/stats/final_stats.txt
```

To run the Nextflow QC and trimming workflow:

```bash
nextflow run nextflow_examples/src/fastpc.nf \
  --fqs_dir /path/to/raw_fastq \
  --qc_out nextflow_examples/results/fastqc \
  --fastp_dir /path/to/trimmed_fastq \
  --cores 4
```
