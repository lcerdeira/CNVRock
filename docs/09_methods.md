# 9. Methods (parameter choices)

This page documents **non-obvious parameter choices** and **why** — the bits
a reviewer will ask about.

## Mapping quality threshold (MQ ≥ 20)

GATK's `--minimum-mapping-quality 20` is the standard threshold for variant
calling (GATK HaplotypeCaller, DELLY, samtools defaults all use it as their
floor or default).

We **lowered from MQ = 40** (the Phase 1 value) because the extended Phase-D
reference contains multiple plasmid contigs that share AMR-cassette sequence
with each other. A read derived from `blaKPC-2` in a sample maps equally well
to several plasmid contigs in the reference, gets MQ = 0, and is discarded
entirely at MQ = 40. Result: `blaKPC-2`, `blaCTX-M-15`, `blaTEM-1` and
`aac6-Ib-cr` returned **zero counts across all 5,000 samples** in the first
exp 32 run.

MQ = 20 corresponds to ≤ 1 % mis-placement probability and recovers the
*near-unique* mappings (those with a discriminating context window). MQ = 10
was considered but rejected as it is rarely used in published pipelines.

We document this clearly in the manuscript Methods section.

## Concurrency cap (10 simultaneous `ascp`)

EBI's Aspera SSH endpoint refuses more than ~10–15 concurrent connections
per source IP. We measured:

| Concurrency | Failure rate |
|---|---|
| 1 | 0 % |
| 10 | 0 % |
| 20 | ≈ 100 % R2 fail |
| 50 | ≈ 100 % fail |

`hpc/aspera_subset_pipeline.sh` is submitted with `--array=…%10` to cap
running tasks at 10. SLURM `MaxArraySize = 5000` is worked around by
chunking the manifest into four ≤ 5000-task arrays with `BATCH_OFFSET`.

## Stratification: species × Bla_Carb × ST cap

See {doc}`04_subset_selection`. ST cap is **150 samples per ST per stratum**
to prevent ST258 / ST11 from dominating the smaller tiers. The full-pool
("80K") tier bypasses the cap.

Sampling within strata is **weighted 1.5× toward carbapenemase carriers**.

## VAE β-warmup

Linear ramp from β = 0 to β = `max_beta` = 1.0 over the first 20 epochs,
then constant. Without β-warmup the KL term dominates early training and
the latent collapses to N(0,1) regardless of input — model produces flat
reconstructions.

## CNV-pattern auxiliary loss

Weight 1.0, warmup 30 epochs. This is a sine-wave-inspired loss that
penalises latent representations failing to track depth gradients across
copy-number boundaries. Empirically required to break the trivial solution
where the VAE encodes only average per-sample depth and reconstructs a flat
profile.

## HMM 6 states with self-transition 0.80

State means initialised at CN ∈ {0, 0.5, 1, 1.5, 2, 3}. Self-transition
0.80 ⇒ expected segment length ≈ 5 kb at 1 kb resolution. Lower
self-transitions over-segment the chromosome; higher values miss small
events (e.g. transposon-bounded gene amplifications).

`hmm_low_cov_threshold = 10` masks bins with < 10 reads — below that the
Gaussian assumption breaks down and the bin is re-imputed from its
neighbours.

## Per-gene PCN thresholds

`pcn_absent_threshold = 0.20`, `pcn_amp_threshold = 1.50` are defaults; the
**per-gene `absent_threshold`** in `plasmid_gene_coords.tsv` overrides
per-gene where the global default would mis-call (e.g. blaCTX-M-15 at 0.50,
blaCTX-M-14 at 1.00 because of higher cross-mapping).

## Chromosomal CRR thresholds

`cnv_crr_amp_threshold = 1.75`, `cnv_crr_gate_threshold = 1.75`. A
gene-region CRR > 1.75 is called as amplified. The gate threshold also
1.75 means a single threshold for both the amp call and the
"signal-strong-enough-to-call" gate. `cnv_min_cn1_proportion = 0.55`
requires the rest of the chromosome to mostly be CN=1 before trusting
the gene-region call (rejects whole-chromosome aneuploidy artefacts).

## Reproducibility seed

The subset selector uses `numpy.random.default_rng(seed=42)`. The same seed
produces identical manifests every time, so the scaling study is fully
deterministic from the same Kleborate ground-truth input.
