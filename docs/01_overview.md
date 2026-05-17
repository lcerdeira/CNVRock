# 1. Overview

## Problem

AMR in *Klebsiella pneumoniae* species complex (KpSC) is driven by two
mechanisms that conventional gene-presence callers (AMRFinder+, Kleborate,
ResFinder) **don't fully capture**:

1. **Chromosomal amplification.** Genes such as `blaSHV` are encoded on the
   KpSC chromosome at single copy in the reference genome, but clinical
   isolates can carry **multiple chromosomal copies** that drive elevated
   resistance. Presence-callers see "blaSHV present" and stop.
2. **Plasmid copy number (PCN).** Plasmid-borne AMR genes (`blaKPC`, `blaNDM`,
   `blaCTX-M`, `blaOXA-48` …) show large per-plasmid PCN variation across
   isolates, again undetectable from presence calls alone.

Both effects compound: a strain may carry blaSHV with 3× chromosomal copies
*and* blaCTX-M-15 on a 5× plasmid. CNVRock recovers both jointly from a
single short-read sequencing experiment.

## Approach

CNVRock is a three-stage pipeline:

1. **1D Conv-VAE** learns a latent representation of 1 kb-binned read-depth
   tensors across the KpSC chromosome (~5,334 bins on NC_016845.1).
2. **Gaussian HMM** segments the per-sample latent reconstructions into
   discrete copy-number states, producing per-bin CN calls.
3. **Gene-level CNV callers** translate per-bin CN segments into per-gene
   copy-number values, separately for chromosomal genes (using CRR = gene /
   flank copy ratio) and plasmid genes (using PCN = gene depth / chromosomal
   median).

## Architecture

```
                                    ┌─────────────────┐
   FASTQ ──BWA──► BAM ──GATK──►  1 kb bin counts   ─┐
                                    └─────────────────┘ │
                                                        ▼
                                              ┌──────────────────┐
                                              │  NPY store       │
                                              │  (n_samples,     │
                                              │   n_bins)        │
                                              └──────────────────┘
                                                        │
                                                        ▼
                                              ┌──────────────────┐
                                              │  1D Conv-VAE     │ ── reconstructions
                                              │  10-dim latent   │
                                              └──────────────────┘
                                                        │
                                                        ▼
                                              ┌──────────────────┐
                                              │  Gaussian HMM    │ ── per-bin CN segments
                                              │  6 states        │
                                              └──────────────────┘
                                                        │
                                          ┌─────────────┴─────────────┐
                                          ▼                           ▼
                                ┌──────────────────┐        ┌──────────────────┐
                                │ Chrom CNV caller │        │ Plasmid CNV      │
                                │ CRR thresholding │        │ PCN thresholding │
                                └──────────────────┘        └──────────────────┘
                                          │                           │
                                          ▼                           ▼
                                ┌──────────────────────────────────────────┐
                                │   evaluation.txt (MCC, FNR, PPV, …)      │
                                └──────────────────────────────────────────┘
```

## Why the scaling study

Phase 1 trained on ~545 KpSC samples. The expansion cohort (this work) makes
~78,000 KpSC samples with FASTQs available. The manuscript reports a **5-point
scaling study** (545 / 5K / 10K / 20K / 40K) showing how chromosomal and
plasmid AMR detection MCC, FNR and PPV change with training-set size.

The hypothesis is **monotonic improvement with diminishing returns** — finding
the saturation point matters both scientifically (what's the cost-effective
training-set size?) and operationally (do we keep scaling to all 78K, or stop
earlier?).
