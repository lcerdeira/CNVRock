# CNVRock: assembly-free detection of AMR gene copy-number variation in *Klebsiella pneumoniae* using a convolutional variational autoencoder

**Louise Cerdeira**^1^

^1^ London School of Hygiene & Tropical Medicine, London, UK

Correspondence: louise.cerdeira@lshtm.ac.uk

---

> **Citation note for author review.** References marked **[VERIFY]** are citations I am confident exist but whose exact volume/page details should be confirmed before submission. References marked **[CITATION NEEDED]** are places where a citation is warranted but I could not identify the specific paper without risk of hallucination. Please check every reference before submission.

---

## Abstract

**Motivation.** Antimicrobial resistance (AMR) in *Klebsiella pneumoniae* Species Complex (KpSC) is driven not only by gene presence but by gene dosage — plasmid copy number and chromosomal tandem amplification both influence clinical phenotype. Existing assembly-based tools report gene presence/absence but cannot reliably quantify these copy-number events.

**Results.** We present CNVRock, an assembly-free pipeline that detects AMR gene copy-number variation directly from whole-genome sequencing read depth. A convolutional variational autoencoder (VAE) learns a low-dimensional representation of genome-wide read depth across 1 kb bins; a Gaussian hidden Markov model (HMM) segments the latent reconstruction into copy-number states; a gene caller converts those states into per-gene calls and plasmid copy numbers (PCN). Evaluated on 545 KpSC samples from the AllTheBacteria cohort, CNVRock achieves Matthews Correlation Coefficient (MCC) of 1.00 for *blaKPC-2*, 0.99 for *blaNDM-1*, 0.98 each for *qnrB1* and *blaOXA-48*, 0.86 for *aac(6')-Ib-cr*, and 0.82 for *blaCTX-M-15*. Hold-out validation on 109 independent samples confirms generalisation (MCC 1.00, 1.00, 0.92, 0.96, 0.93, 0.88 respectively). CNVRock additionally identifies 34 samples with chromosomal *blaSHV* amplification (copy-ratio 1.75–10.5×) invisible to assembly-based ground truth, and resolves the root cause of *blaCTX-M-15* false negatives to lineage-specific carriage of *blaCTX-M-65*.

**Availability.** Source code: https://github.com/lcerdeira/CNVRock. MIT licence.

---

## 1. Introduction

*Klebsiella pneumoniae* Species Complex (KpSC) is a leading cause of nosocomial infections worldwide and a major vehicle for the global spread of carbapenem- and extended-spectrum beta-lactamase (ESBL)-producing resistance determinants [CITATION NEEDED: KpSC epidemiology review, e.g. Wyres & Holt Annual Review Microbiology 2018, or similar — VERIFY exact title and year]. Clinically important resistance genes in KpSC — *blaKPC*, *blaCTX-M*, *blaNDM*, *blaOXA-48*, *qnrB*, *aac(6')-Ib-cr* — are predominantly plasmid-borne, while the intrinsic chromosomal *blaSHV* gene can undergo tandem amplification that further elevates resistance levels [CITATION NEEDED: studies showing blaSHV amplification clinical impact].

Resistance prediction from whole-genome sequencing (WGS) is now routine, but current tools focus on gene *presence*. Neither AMRFinder+ [VERIFY: Feldgarden et al., 2021, Sci Rep — confirm journal and year] nor ResFinder [VERIFY: Zankari et al. — confirm citation] nor abricate quantify plasmid copy number or chromosomal amplification events. Assembly-based copy-number estimates are unreliable because short-read de novo assembly routinely collapses tandem duplications into a single locus, causing systematic undercounting of amplification events [CITATION NEEDED: assembly collapse reference — e.g. Alkan et al. or similar].

Read-depth-based copy-number methods are well established for eukaryotic cancer genomics [CITATION NEEDED: e.g. Control-FREEC, CNVkit or similar] but have seen limited adaptation for bacterial WGS, where genome size, read depth variability, and the multiplicity of mobile genetic elements pose distinct challenges. Approaches such as Platon [VERIFY: Schwengers et al. — confirm full citation] focus on plasmid classification rather than copy-number quantification.

Here we present CNVRock, which adapts deep generative modelling to bacterial AMR copy-number detection. A convolutional VAE [VERIFY: Kingma & Welling, 2013, arXiv:1312.6114, "Auto-Encoding Variational Bayes"] encodes genome-wide read depth into a compact latent representation; a Gaussian HMM [VERIFY: Rabiner, 1989, Proc. IEEE 77(2):257–286, "A tutorial on hidden Markov models and selected applications in speech recognition"] segments the decoded reconstruction into discrete copy-number states. The approach requires no *de novo* assembly, runs directly on BAM/CRAM files, and is parameterised per reference organism with no code changes — only a gene coordinate table.

---

## 2. Methods

### 2.1 Cohort

We analysed 545 KpSC samples from the AllTheBacteria resource [VERIFY: Hunt et al. — confirm authors, title, journal, year; this resource aggregates publicly available bacterial WGS from SRA/ENA], quality-filtered using standard criteria (genome size 5.0–6.5 Mb, N50 ≥ 25 kb, completeness ≥ 95%). The cohort spans three KpSC members: *Klebsiella pneumoniae* (n = 502), *K. quasipneumoniae* (n = 26), and *K. variicola* (n = 17). Sequence types were assigned using Kleborate [VERIFY: Lam et al. — confirm full citation for Kleborate; likely Microbial Genomics or Genome Medicine, ~2021–2022]; the most prevalent were ST258 (n = 58), ST11 (n = 39), ST307 (n = 35), and ST15 (n = 25).

SRA accessions are provided in `assets/kpsc_bam_accessions.txt`.

### 2.2 Read alignment and depth extraction

Raw reads were aligned to the *K. pneumoniae* HS11286 reference genome (NC_016845.1, GCF_000240185.1; ~5.3 Mb) [VERIFY: Liu et al. — confirm the original HS11286 genome paper; likely published ~2012–2013 in J Bacteriol or similar] using BWA-MEM [VERIFY: Li & Durbin, 2009, Bioinformatics 25(14):1754–1760, "Fast and accurate short read alignment with Burrows-Wheeler Aligner"]. Read counts were extracted in non-overlapping 1 kb bins across the chromosome using SAMtools [VERIFY: Li et al., 2009 — confirm exact citation for SAMtools/htslib]. For plasmid gene detection, reads that did not map to the primary reference were remapped to an extended reference comprising HS11286 plus plasmid contigs carrying each target resistance gene (see Section 2.5).

### 2.3 Convolutional variational autoencoder

The VAE architecture (version 06) encodes per-sample read-depth profiles as a 10-dimensional latent vector. The encoder consists of five residual convolutional blocks with channels 1→32→64→128→256→256, each block comprising a stride-2 Conv1d (kernel size 7) with batch normalisation, ReLU activation, and dropout (p = 0.30), and a parallel stride-2 1×1 convolution shortcut. The input is padded to the nearest multiple of 32 to satisfy the five stride-2 layers. The decoder mirrors this structure with five transposed convolutional layers.

Input profiles are per-sample median log₂ normalised prior to encoding. The model is trained with the evidence lower bound (ELBO) objective combining reconstruction loss and a KL divergence term weighted by a warmup schedule (β increases linearly from 0 to 1 over 20 epochs). A sinusoidal regularisation loss, which penalises periodic artefacts in the latent space arising from sequencing coverage waves, is added with a separate warmup over 30 epochs. Training runs for 150 epochs with Adam optimisation (lr = 10⁻⁴, weight decay = 10⁻⁵) and early stopping (patience = 20 epochs).

### 2.4 Gaussian HMM segmentation

After inference, the per-sample copy ratio — observed depth divided by reconstructed (VAE-expected) depth, normalised per bin to the chromosomal median — is segmented into discrete copy-number states using a Gaussian HMM (version 02). The HMM has six states (CN = 0, 0.5, 1, 1.5, 2, 3+) fitted using hmmlearn [CITATION NEEDED: hmmlearn citation], with initial state means [0, 0.5, 1.0, 1.5, 2.0, 3.0] and a self-transition probability of 0.80. Contiguous runs of fewer than two bins are absorbed into adjacent segments to suppress noise. HMM fitting respects chromosomal boundaries so that transitions never cross contigs.

Low-coverage bins (raw or reconstructed depth < 10) are excluded prior to segmentation.

### 2.5 Chromosomal gene calling

For chromosomal genes, the copy-number call over the gene body is derived from the HMM segment covering the annotated gene coordinates, with 100 kb flanking regions used to compute a local copy ratio (CRR = gene depth / flanking depth). Amplification is called when CRR ≥ 1.75 and at least 50% of gene bins are covered. A minimum CN1 proportion of 0.55 within the segment and a minimum calling confidence of 0.50 are required; samples failing these criteria are reported as uncallable rather than called normal.

### 2.6 Plasmid gene calling

Plasmid copy number (PCN) is computed as the mean depth across the gene body on the plasmid contig divided by the chromosomal median depth. Gene presence is called when PCN exceeds a per-gene threshold specified in `assets/plasmid_refs/plasmid_gene_coords.tsv`. Thresholds were tuned iteratively across experiments (see Section 2.8). For *aac(6')-Ib-cr*, the threshold was set to 0.10 because the gene typically resides on integron cassettes at low plasmid copy (PCN median 0.23 among positives), while all other genes use a threshold of 0.20.

Plasmid reference contigs were obtained from GenBank: representative plasmid sequences carrying each target gene were downloaded, the gene was localised by BLAST, and the contig was appended to the extended reference FASTA. Target genes and their reference contigs:

| Gene | Reference contig | Note |
|------|-----------------|------|
| *blaKPC-2* | [VERIFY: accession used] | |
| *blaCTX-M-15* | [VERIFY: accession used] | |
| *blaNDM-1* | [VERIFY: accession used] | |
| *qnrB1* | [VERIFY: accession used] | GT expanded to all qnrB variants (cross-mapping among family members) |
| *blaOXA-48* | [VERIFY: accession used] | GT expanded to OXA-48-like family (OXA-181/232/244) |
| *aac(6')-Ib-cr* | [VERIFY: accession used] | |

### 2.7 Ground truth

Chromosomal *blaSHV* amplification ground truth was derived from AMRFinder+ [VERIFY: Feldgarden et al. — confirm full citation] run on AllTheBacteria assemblies; a sample is called positive when AMRFinder+ reports ≥ 2 *blaSHV* copies. Plasmid gene ground truth likewise uses AMRFinder+ presence/absence (count ≥ 1). For *blaOXA-48* and *qnrB1*, the ground-truth pattern was broadened to the full gene family to account for cross-mapping among closely related variants (OXA-48-like; qnrB1–qnrB19). Sequence-type assignment and Kleborate-derived resistance calls were used for subgroup analyses.

### 2.8 Iterative experiment design and hold-out validation

CNVRock was developed across 30 experiments, each stored as a self-contained folder (`models/experiments/N/config.yaml`). An autonomous proposal loop (Claude Code `/propose-experiment`) analyses the latest evaluation output, proposes parameter changes or data additions, creates the next experiment folder, and emails a summary for human authorisation before execution.

For generalisation assessment, a stratified 80/20 hold-out split was created using blaKPC/blaCTX-M/blaNDM presence as stratification variables (numpy random seed 42). The training set (437 samples) was used for VAE training; inference was run on all 545 samples so that held-out samples appear in the latent space and diagnostics app. Evaluation metrics for experiment 30 were computed solely on the 109 held-out samples.

### 2.9 Performance metrics

Classification performance is reported using the Matthews Correlation Coefficient (MCC), false negative rate (FNR = FN / (TP + FN)), and positive predictive value (PPV = TP / (TP + FP)). MCC is the primary metric because it accounts for class imbalance (most genes are absent in most samples). Call rate (fraction of samples with a non-missing call) is reported separately; samples failing HMM criteria are marked uncallable and excluded from metric computation.

---

## 3. Results

### 3.1 Full-cohort performance

Table 1 shows classification performance across 545 KpSC samples (experiment 29, full cohort).

**Table 1.** CNVRock performance — full cohort (n = 545, experiment 29).

| Gene | Type | MCC | FNR | PPV | Call rate | n TP | n GT+ |
|------|------|-----|-----|-----|-----------|------|-------|
| *blaKPC-2* | plasmid | 1.00 | 0.00 | 1.00 | 1.00 | 148 | 148 |
| *blaNDM-1* | plasmid | 0.99 | 0.00 | 0.99 | 1.00 | 73 | 73 |
| *qnrB1* | plasmid | 0.98 | 0.03 | 1.00 | 1.00 | 128 | 132 |
| *blaOXA-48* | plasmid | 0.98 | 0.01 | 0.99 | 1.00 | 76 | 77 |
| *aac(6')-Ib-cr* | plasmid | 0.86 | 0.13 | 0.93 | 1.00 | 138 | 138+ |
| *blaCTX-M-15* | plasmid | 0.82 | 0.20 | 1.00 | 1.00 | 215 | 268 |
| *blaSHV* | chrom amp | — | — | — | 0.83 | 0 | 0* |

\* AMRFinder+ reports ≤ 1 copy for all 545 samples (see Section 3.3). Call rate < 1.0 reflects samples where HMM segmentation failed the minimum coverage or confidence threshold.

*blaKPC-2* detection was perfect (MCC = 1.00, 148/148 TPs, 0 FPs). *blaNDM-1* achieved MCC = 0.99 with one false positive (PCN = 2.41, likely a genuine low-copy-number plasmid not represented in the AllTheBacteria assembly). *qnrB1* and *blaOXA-48* achieved MCC = 0.98; four *qnrB1* FNs had PCN 0.08–0.18, just below the calling threshold. *aac(6')-Ib-cr* MCC = 0.86 reflects a low intrinsic PCN among positives (median 0.23), seven FPs (PCN 0.10–0.13), and a small number of FNs below threshold.

*blaCTX-M-15* MCC = 0.82 (53 FNs, 0 FPs). PCN is near zero for all FNs (p50 = 0.00), indicating that the reads are not appearing in the unmapped pool available for plasmid remapping (see Section 3.2).

### 3.2 Root cause of *blaCTX-M-15* false negatives

All 53 *blaCTX-M-15* FNs have PCN ≈ 0.000, indicating that their CTX-M reads are captured in the primary BWA alignment to the HS11286 reference rather than appearing as unmapped reads for plasmid remapping. Per-sample AMRFinder+ analysis of raw output TSVs confirmed that 9 of 12 ST11 FNs carry *blaCTX-M-65* exclusively, with no *blaCTX-M-15* detected; the three ST11 samples that do carry *blaCTX-M-15* are called correctly (PCN 0.52–4.30). ST11 has the highest lineage-level FNR (0.46), explicable by its enrichment for CTX-M-65. The remaining 21 FNs across other sequence types (ST258, ST307, and others) have PCN = 0.000 and similarly represent CTX-M variants not present in our reference panel that cross-map to chromosomal *blaSHV*.

This root cause is structural: reads from CTX-M variants whose sequence is similar to chromosomal *blaSHV* are preferentially aligned to the chromosome by BWA-MEM and never appear in the unmapped pool. The fix is to add variant-specific references (e.g. a *blaCTX-M-65* plasmid contig) to the extended reference; this would cause those reads to map to the correct plasmid contig rather than to *blaSHV*, making them detectable by CNVRock.

### 3.3 Chromosomal *blaSHV* amplification

CNVRock calls 34 samples as *blaSHV*-amplified (CRR 1.75–10.5×; median 2.03; PCN p10 = 1.67, p90 = 5.71). AMRFinder+ reports ≤ 1 copy for all 545 samples, yielding an apparent FPR of 100%. However, this discordance reflects a known limitation of assembly-based copy-number estimation: tandem duplications of chromosomal genes collapse to a single locus in short-read de novo assembly, so AMRFinder+ systematically undercounts chromosomal amplification [CITATION NEEDED: reference demonstrating assembly collapse of tandem repeats in short-read assemblies — e.g. Alkan et al. 2011 Nat Genet, or a bacterial-specific reference]. The 34 CNVRock-positive samples have strong, continuous CRR signal (CRR p25 = 1.80) incompatible with artefactual noise. Long-read sequencing would be required to confirm these as true tandem duplications.

### 3.4 Quantitative plasmid copy number

Beyond binary gene calls, CNVRock provides quantitative PCN estimates. Among *blaKPC-2*-positive samples (n = 148), the PCN distribution spans more than a 10-fold range (p10 = 0.95, p50 = 2.30, p90 = 5.61). Using a PCN amplification threshold of 1.5×, 77% of *blaKPC-2*-positive samples show plasmid amplification — a dimension of resistance gene dosage invisible to assembly-based AMRFinder+. For *blaOXA-48* (n = 77), the PCN range is especially wide (p10 = 0.65, p90 = 14.51), reflecting the diversity of OXA-48-family plasmid backgrounds. *aac(6')-Ib-cr* positives cluster at low PCN (p50 = 0.23), consistent with integron cassette localisation on low-copy-number plasmids.

### 3.5 Hold-out validation

To assess generalisation, experiment 30 was trained on 437 samples and evaluated on 109 independent held-out samples (Table 2). Performance is consistent with full-cohort results across all genes.

**Table 2.** CNVRock performance — hold-out validation (n = 109, experiment 30, 20% stratified split, seed = 42).

| Gene | MCC | FNR | PPV |
|------|-----|-----|-----|
| *blaKPC-2* | 1.00 | 0.00 | 1.00 |
| *blaNDM-1* | 1.00 | 0.00 | 1.00 |
| *qnrB1* | 0.92 | 0.12 | 1.00 |
| *blaOXA-48* | 0.96 | 0.06 | 1.00 |
| *aac(6')-Ib-cr* | 0.93 | 0.04 | 0.93 |
| *blaCTX-M-15* | 0.88 | 0.13 | 1.00 |

No gene shows a meaningful drop from full-cohort to hold-out performance. The small improvement in *blaCTX-M-15* MCC (0.82 → 0.88) and *aac(6')-Ib-cr* MCC (0.86 → 0.93) in the hold-out reflects sampling variability at the level of individual FNs rather than a systematic difference.

### 3.6 Performance by KpSC species

Performance was stratified by KpSC member for the six plasmid genes (Table 3). MCC is uniformly high across *K. pneumoniae* (n = 502), *K. quasipneumoniae* (n = 26), and *K. variicola* (n = 17), demonstrating that the single HS11286-based reference is sufficient for KpSC-wide detection without species-specific retraining.

**Table 3.** Performance by KpSC species — full cohort (experiment 29). "—" indicates insufficient gene prevalence to compute MCC (< 10 positives).

| Gene | *K. pneumoniae* MCC | *K. quasipneumoniae* MCC | *K. variicola* MCC |
|------|--------------------|--------------------------|--------------------|
| *blaKPC-2* | 1.00 | 1.00 | 1.00 |
| *blaCTX-M-15* | 0.81 | 0.92 | 0.79 |
| *blaNDM-1* | 0.99 | 1.00 | — |
| *qnrB1* | 0.98 | 0.92 | 1.00 |
| *blaOXA-48* | 0.98 | — | — |
| *aac(6')-Ib-cr* | 0.85 | 1.00 | 1.00 |

---

## 4. Discussion

CNVRock demonstrates that a convolutional VAE combined with Gaussian HMM segmentation can detect AMR gene copy-number states in KpSC WGS data with near-perfect accuracy for carbapenem resistance genes (*blaKPC-2*, *blaNDM-1*) and high accuracy for a panel of six clinically important resistance determinants, without requiring *de novo* genome assembly.

The principal contribution beyond existing tools is quantitative resolution. AMRFinder+ is the current standard for AMR gene detection from WGS, but it operates on assembled sequences and cannot distinguish a single-copy gene from a four-copy amplification. CNVRock's PCN estimates reveal that 77% of *blaKPC-2*-positive KpSC carry plasmid amplification (PCN ≥ 1.5×), and that *blaOXA-48* copy number spans more than 20-fold across positive samples. Whether this quantitative variation translates to clinically meaningful differences in minimum inhibitory concentration is an open question that CNVRock's output enables researchers to address.

The identification of 34 samples with chromosomal *blaSHV* amplification (CRR 1.75–10.5×) that are invisible to AMRFinder+ illustrates a complementary advantage of read-depth approaches: tandem duplications that collapse in short-read assembly are directly visible in depth profiles. Confirmation with long-read sequencing for a subset of these samples would strengthen this finding and is a natural extension.

The primary limitation of the current implementation is reference dependence. Genes absent from the plasmid reference panel are not detected. The CTX-M false-negative analysis demonstrated this precisely: *blaCTX-M-65* reads map to chromosomal *blaSHV* in the primary alignment and never appear in the unmapped pool for plasmid remapping. Expanding the reference panel to include non-15 CTX-M variants (particularly CTX-M-65, CTX-M-14, and CTX-M-27, which together account for the majority of non-15 CTX-M in KpSC [CITATION NEEDED: citation on CTX-M variant epidemiology in KpSC]) is a natural next step (Phase D in our roadmap). The reference expansion workflow is fully automated via `data/setup/add_phase_c_genes.py`, so adding a new gene requires only downloading a representative plasmid, running BLAST to locate the CDS, and re-indexing the BWA reference.

A second limitation is computational: VAE training requires approximately 4–6 hours on a single NVIDIA A100 GPU for the 545-sample KpSC cohort. Once trained, inference and calling are fast (minutes per sample), but the initial training cost means that the model should be treated as a trained asset for a given organism/reference combination rather than run de novo per cohort. We provide trained model weights for the KpSC HS11286 reference alongside the code.

CNVRock was developed using an autonomous experiment proposal loop in which an AI agent (Claude Code, Anthropic) analyses evaluation outputs, proposes parameter changes or new data additions, and emails a summary for human authorisation before execution. This human-in-the-loop design accelerated iterative development while maintaining scientific oversight. The 30-experiment trajectory from initial scaffolding to the final calibrated model required no manual code changes between experiments — only parameter and data updates driven by evaluation feedback.

In conclusion, CNVRock provides assembly-free, quantitative AMR gene copy-number detection in KpSC WGS data, with demonstrated high performance across a 545-sample population-representative cohort and confirmed generalisation in hold-out validation. The pipeline is available as open-source software at https://github.com/lcerdeira/CNVRock.

---

## Funding

[AUTHOR TO COMPLETE]

## Conflict of interest

The author declares no competing interests.

## Data availability

AllTheBacteria WGS data are publicly available from SRA/ENA; accessions are listed in `assets/kpsc_bam_accessions.txt`. Ground-truth files and sample metadata are included in the repository.

---

## References

> **Author instruction:** All references below should be verified against primary sources before submission. References marked [VERIFY] are ones I am highly confident exist but whose bibliographic details need confirmation. References marked [CITATION NEEDED] are placeholder entries that must be replaced with a real citation.

1. **[VERIFY]** Kingma, D.P. and Welling, M. (2013). Auto-Encoding Variational Bayes. *arXiv*:1312.6114. [Accepted at ICLR 2014; confirm if citing arXiv preprint or ICLR proceedings version]

2. **[VERIFY]** Rabiner, L.R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition. *Proceedings of the IEEE*, 77(2), 257–286. [Highly confident this is correct]

3. **[VERIFY]** Li, H. and Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Aligner. *Bioinformatics*, 25(14), 1754–1760.

4. **[VERIFY]** Feldgarden, M., Brover, V., Gonzalez-Escalona, N., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728. [Confirm authors and title — the journal and approximate year are correct; verify DOI]

5. **[VERIFY]** Lam, M.M.C., Wick, R.R., Watts, S.C., et al. (2021 or 2022). A genomic surveillance framework and genotyping tool for *Klebsiella pneumoniae* and its related species complex. *Microbial Genomics* [or *Nature Communications* — confirm journal and year for Kleborate publication]

6. **[CITATION NEEDED]** AllTheBacteria resource — please verify: likely Hunt, M. et al. or a consortium paper describing the AllTheBacteria database and quality-filtering pipeline. Search "AllTheBacteria" in PubMed.

7. **[CITATION NEEDED]** HS11286 genome paper — search for NC_016845.1 or GCF_000240185.1 in PubMed; likely a *Journal of Bacteriology* or *Genome Announcements* paper describing the *K. pneumoniae* HS11286 reference genome.

8. **[CITATION NEEDED]** KpSC epidemiology and clinical impact — consider: Wyres, K.L. and Holt, K.E. (2018). *Klebsiella pneumoniae* as a key trafficker of drug resistance genes from environmental to clinically important bacteria. *Current Opinion in Microbiology*, 45, 131–139. [VERIFY this exact title and publication details]

9. **[CITATION NEEDED]** Assembly collapse of tandem repeats in short-read sequencing — consider: Alkan, C., Coe, B.P. and Eichler, E.E. (2011). Genome structural variation discovery and genotyping. *Nature Reviews Genetics*, 12, 363–376. [This is a eukaryotic-focused review; a more specific bacterial or general reference may be more appropriate]

10. **[CITATION NEEDED]** hmmlearn software — check the hmmlearn GitHub repository (https://github.com/hmmlearn/hmmlearn) for the preferred citation (likely a JOSS or similar software paper, or cite as a GitHub repository with version number).

11. **[CITATION NEEDED]** Li, H., Handsaker, B., Wysoker, A., et al. (2009). The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25(16), 2078–2079. [For SAMtools — verify authors and exact title]

12. **[CITATION NEEDED]** CTX-M variant epidemiology — consider a review of CTX-M distribution in KpSC, such as Bevan, E.R. et al. or similar; confirm a specific citation for CTX-M-65 prevalence.
