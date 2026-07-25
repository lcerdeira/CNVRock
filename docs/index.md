# CNVRock

**AMR gene copy-number variation in *Klebsiella pneumoniae* species complex (KpSC) using variational autoencoders.**

Antimicrobial resistance (AMR) in KpSC is driven not just by gene **presence**
but by gene **copy number** — amplification of chromosomal AMR loci and variable
plasmid copy number (PCN) both have direct clinical consequences. CNVRock
combines a 1D convolutional VAE, a Gaussian HMM segmenter, and gene-level CNV
callers to recover both chromosomal and plasmid CN events from short-read
sequencing.

This site documents the **end-to-end pipeline** used in the manuscript:
sample selection, data acquisition, reference preparation, model training, the
scaling study (5K → 10K → 20K; the pre-registered 40K and 80K tiers were
deprecated once performance was shown to saturate by 10K), and the
[validation and robustness](11_validation.md) analyses that test the
copy-number values themselves.

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

09_methods.md
10_reproducibility.md
11_validation.md
```

## Quick links

- **Repository:** <https://github.com/lcerdeira/CNVRock>
- **Interactive demo:** <https://cnvrock.streamlit.app/> — sample-viewer (per-sample copy-number profile + latents), Monitor, and latent-coverage pages, on bundled 200-sample subsamples of all three organisms (*K. pneumoniae* 10 K tier, *A. baumannii*, *C. auris*). No install required. Each dataset surfaces **showcase isolates** with the strongest signals (e.g. *C. auris* ERG11 15× amplification, chr5 aneuploidy; *K. pneumoniae* blaSHV 85× tandem; *A. baumannii* blaOXA-23 80×). The bundles embed the read-count inputs store so the copy-number viewer runs on the cloud; rebuild them with `diagnostics/build_demo_bundle.py --store … --showcase …`.
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
`main`. The KpSC scaling study and the evolutionary analyses (selection
dynamics, genome-wide scan, co-occurrence) now live under `kpsc_paper/`,
where they are being prepared as a separate manuscript.
```
