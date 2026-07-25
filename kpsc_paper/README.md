# KpSC copy-number evolution — second manuscript

Material split out of the main CNVRock methods paper on 2026-07-22, after
reviewers judged the *Klebsiella pneumoniae* species complex (KpSC) analyses
to read as a population-genetics/evolution study rather than an AMR
copy-number methods study, and asked for them to be removed.

That judgement is defensible: the analyses collected here ask evolutionary
questions (is amplification under selection? does it sweep? does it track
sequence type and calendar time?), whereas the methods paper asks a
measurement question (can we quantify AMR gene dosage from short reads,
across organisms?). They are two papers.

## What moved here

| Script | Question |
|---|---|
| `analysis/birth_death_cnv.py` | lambda/mu birth-death estimates — is amplification under selection? |
| `analysis/genome_wide_cnv_scan.py` | genome-wide copy-ratio x phenotype scan |
| `analysis/annotate_scan_candidates.py`, `scan_candidate_mge.py` | functional / mobile-element annotation of scan hits |
| `analysis/cnv_cooccurrence_network.py`, `cnv_cooccurrence_v2.py` | plasmid exclusion and co-amplification modules |
| `analysis/cnv_st_mixedmodel.py` | CNV x sequence-type association |
| `analysis/temporal_amplification.py` | calendar-time trend in amplification prevalence |
| `analysis/kvariicola_multiref.py` | K. variicola blaLEN multi-reference calling |
| `analysis/porin_cnv_analysis.py` | KpSC porin copy number |
| `hpc/*` | scan, Kleborate GT and K. variicola HPC jobs |
| `docs/08_scaling_study.md` | 5K/10K/20K scaling study |

## What deliberately did NOT move

The following were run on KpSC data but are **methodological validation**, not
KpSC biology, and stay with the methods paper:

- `analysis/longread_depth_validation.py` — the ONT long-read agreement
  (Pearson r = 0.96, n = 255). This is the only external validation of the
  copy-number *values* anywhere in the project. Removing it would leave the
  methods paper with no orthogonal validation at all.
- `analysis/vae_ablation.py` — the three-way baseline ablation (genome-wide
  median vs 200 housekeeping bins vs Conv-VAE). It exists because an earlier
  reviewer asked whether the learned baseline earns its place, and it is the
  direct evidence bearing on the current reviewer's proposal to normalise on
  7 MLST loci. The housekeeping baseline lost on all three metrics
  (noise floor 0.151 vs 0.112; spike-in RMSE 1.038 vs 0.710; false-amp rate
  0.031 vs 0.010) — and 7 loci is a sparser version of that baseline.
- `analysis/percall_uncertainty.py`, `threshold_sensitivity.py`,
  `mc_dropout_uncertainty.py`, `latent_intrinsic_dim.py` — per-call CIs,
  threshold calibration, uncertainty and latent-dimension ablation.

The honest framing for the methods paper is that the pipeline was **developed
and validated** on a KpSC cohort and then **applied** to four other organisms.
Reviewers object to a Klebsiella *narrative*, not to Klebsiella existing as a
development set.

## Entanglement to resolve

The E. coli extended reference reuses 17 KpSC plasmid contigs (only
CP061207.1 is E. coli-native). Four genes — blaTEM, sul1, sul2, qnrB — lose
sensitivity because E. coli variant reads map inefficiently onto Klebsiella
backbones. If the methods paper must be free of KpSC assets, those contigs
need E. coli-native replacements, which would also fix the weak detection.

## Status

Parked. The scaling study, scan, selection dynamics and co-occurrence results
are complete and reproducible; the writing is not started. Text for these
analyses currently lives in the main manuscript and needs extracting.
