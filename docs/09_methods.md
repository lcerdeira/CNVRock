# 9. Methods (parameter choices)

This page documents **non-obvious parameter choices** and **why** — the bits
a reviewer will ask about.

## Mapping quality threshold (MQ ≥ 0) and gene-family aggregation

This is the parameter choice that took the longest to get right. The
diagnostic narrative is below; the final settings are:

- GATK `CollectReadCounts --minimum-mapping-quality 0`
- Plasmid PCN evaluated by **gene family** (sum of bin counts across all
  allele-variant CDSs in each family), using
  `assets/plasmid_refs/plasmid_gene_families.tsv`.

### Why MQ = 0

We first ran the pipeline at **MQ = 40** (the Phase 1 value), which produced
**zero plasmid-gene counts** for `blaKPC-2`, `blaCTX-M-15`, `blaTEM-1`,
`aac6-Ib-cr` across all 5,000 samples — while `qnrB1`, `blaOXA-181`,
`blaNDM-5` were fine. We assumed multi-mapping with too-strict MQ and
relaxed to MQ = 20 (the GATK / DELLY / samtools default). It did not help.

Cross-checking the gene coordinates against the curated GenBank annotations
of every contig in `HS11286_extended.fasta` revealed that
`plasmid_gene_coords.tsv` **was correct all along** — the coordinates match
GenBank to within ≤ 3 bp for every gene.

The real cause is **multi-mapping between near-identical allele variants**:

| Family | Member CDSs (in our ref) | nt identity |
|---|---|---|
| `blaNDM` | NDM-1 (MZ606384.2), NDM-5 (CP034201.2) | ~99.8 % (2 aa / 6 nt difference) |
| `blaOXA-48-like` | OXA-48 (JN626286.1), OXA-181 (CP113224.1) | ~98 % (4 aa difference) |
| `blaCTX-M` | CTX-M-15 (MK552109.1), CTX-M-14 (NZ_CP031850.1), CTX-M-65 (CP030319.1), CTX-M-27 (CP138680.1) | ~95–99 % within each subgroup |

A read derived from any one of these variants maps **equally well** to all of
them in the same family. BWA assigns MQ = 0 to these multi-mappers and
places the read at one location at random. At any MQ ≥ 1 they are filtered
out entirely.

**MQ = 0 keeps them.** Because each multi-mapped read still ends up at *one*
location (chosen by BWA), individual variant CDSs receive a fraction of the
true coverage. Summing bin counts across all family members recovers the
full per-sample per-family signal.

### Why this is methodologically sound

1. The multi-mapping is **between annotated allele variants of the same gene
   family** — not between unrelated genomic regions. Counting these reads as
   "the family is present" is biologically correct.
2. **Clinical reporting works at family level.** AMRFinder+ output is
   variant-level but Kleborate-style summary tables report `Bla_Carb_acquired
   = NDM-1` or `NDM-5` interchangeably; treatment recommendations are
   family-level.
3. For chromosomal CNV detection (the `blaSHV` extra-copy story) the KpSC
   chromosome is largely uniquely-mappable. MQ = 0 vs MQ = 20 affects only
   repeat / transposon regions which the VAE absorbs as noise. Single
   pipeline, single count file per sample.
4. The cost of getting MQ = 0 wrong (over-counts from spurious mappings) is
   small: BWA's secondary-alignment quality scores still exclude very poor
   matches, and the 1 kb bin is large enough that occasional spurious
   placements don't dominate the per-sample family count.

### Why gene-family aggregation, not individual genes

Individual blaNDM-1 vs blaNDM-5 detection is fundamentally impossible from
short reads when both reference CDSs are in the alignment target: BWA
*cannot* distinguish them above MQ = 0. AMRFinder+ resolves this differently
— it works on **assemblies**, not raw reads, and exploits *the rest of the
plasmid backbone* as discriminating context. CNVRock works on read depth,
so we accept the limitation and report at family level. The variant call
(if needed) comes from the Kleborate column in the ground-truth join.

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
