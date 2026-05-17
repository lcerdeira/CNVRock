# CNVRock

**AMR gene copy-number variation in *Klebsiella pneumoniae* species complex (KpSC) using variational autoencoders.**

Antimicrobial resistance (AMR) in KpSC is driven not just by gene **presence**
but by gene **copy number** — amplification of chromosomal AMR loci and variable
plasmid copy number (PCN) both have direct clinical consequences. CNVRock
combines a 1D convolutional VAE, a Gaussian HMM segmenter, and gene-level CNV
callers to recover both chromosomal and plasmid CN events from short-read
sequencing.

This site documents the **end-to-end pipeline** used in the manuscript:
sample selection, data acquisition, reference preparation, model training, and
the scaling study (5K → 10K → 20K → 40K).

```{toctree}
:maxdepth: 2
:caption: Pipeline

01_overview.md
02_data_acquisition.md
03_reference_and_intervals.md
04_subset_selection.md
05_npy_stores.md
06_training.md
07_evaluation.md
```

```{toctree}
:maxdepth: 2
:caption: Manuscript

08_scaling_study.md
09_methods.md
10_reproducibility.md
```

## Quick links

- **Repository:** <https://github.com/lcerdeira/CNVRock>
- **Interactive diagnostics (demo):** <https://cnvrock.streamlit.app/> — Monitor / sample viewer / latent-coverage pages on a bundled 200-sample subsample of exp 32. No install required.
- **Manuscript reference cohort:** Phase 2 expansion (~78,000 KpSC samples with FASTQ URLs on ENA)
- **Compute:** LSHTM HPC, GPU partition (A40), Aspera SDK for high-throughput EBI downloads

## At a glance

| Stage | Tool / Asset |
|---|---|
| Stratified sample selection | `data/setup/select_expansion_subset.py` |
| FASTQ download | IBM Aspera (`ascp` via bioconda `aspera-cli`) |
| Alignment | BWA-MEM 0.7.18, paired-end |
| Read counts | GATK `CollectReadCounts`, 1 kb bins, MQ ≥ 20 |
| Reference | `HS11286_extended.fasta` (chrom + 11 plasmid contigs) |
| Training | 1D Conv-VAE + Gaussian HMM, `models/train.py` |
| Plasmid CNV calls | `models/cnv/07_plasmid_cnv_caller.py` |
| Chromosomal CNV calls | `models/cnv/06_gene_cnv_caller.py` |
| Evaluation | AMRFinder+ and Kleborate ground truth |

```{note}
This documentation is built automatically by Read the Docs on every push to
`main`. The scaling-study tables and figures populate as the experiments
complete (see {doc}`08_scaling_study`).
```
