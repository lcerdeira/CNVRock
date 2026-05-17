# 4. Subset selection

For the scaling study we generate **four strictly-nested manifests**:

```
5,000  ⊂  10,000  ⊂  20,000  ⊂  40,000  ⊂  77,906 (full eligible)
```

A sample in a smaller subset is **always also in every larger one**. This is
critical for the scaling study: any change in metrics is attributable to
*more data*, not to a different sample mix.

Script: `data/setup/select_expansion_subset.py` (run once locally — outputs
all four manifests to `assets/`).

## Eligibility pipeline

Starting from 88,128 KpSC samples in our Kleborate ground truth:

```
88,128  Kleborate-typed samples
   ↓    KpSC core species filter
85,339  K. pneumoniae / K. quasipneumoniae / K. variicola / K. africana
   ↓    BioSample → run-accession bridge (kpsc_expansion_metadata.tsv)
85,339
   ↓    inner join with ENA URL manifest (must have FASTQs)
77,906  samples eligible for download
```

## Stratification

Inside the eligible pool we **stratify by species × Bla_Carb_acquired**
(carbapenemase carrier yes/no), and within each stratum **cap each ST at 150
samples** to avoid over-weighting common clones (ST258, ST11, …). This caps
the pool at 38,192 ST-controlled samples; for the **40K and 80K tiers** we
draw additional samples without the ST cap to reach the target size.

Inside each stratum samples are weighted **1.5× toward carbapenemase
carriers** before drawing — they're the positive class CNVRock cares about,
so we want them well-represented in the smaller tiers.

## Nesting trick

We do **one** weighted shuffle of the full eligible pool with `seed=42`. The
manifests are then **prefixes** of that single ordering:

```python
ordered = rng.choice(eligible, size=len(eligible), replace=False, p=weights)
manifest_5k  = ordered[:5_000]
manifest_10k = ordered[:10_000]
manifest_20k = ordered[:20_000]
manifest_40k = ordered[:40_000]
manifest_80k = ordered          # full pool
```

By construction, **every smaller manifest is a strict subset** of every larger
one. Verified at write time:

```bash
$ wc -l assets/kpsc_expansion_subset_*.tsv
   5001   kpsc_expansion_subset_5k.tsv
  10001   kpsc_expansion_subset_10k.tsv
  20001   kpsc_expansion_subset_20k.tsv
  40001   kpsc_expansion_subset_40k.tsv
  77907   kpsc_expansion_subset_80k.tsv
```

## Anchoring to in-flight data

If `assets/kpsc_expansion_subset_5k.tsv` **already exists** when the script is
re-run (e.g. because we already started downloading those samples), the script
**locks in** that 5K and only varies the *remainder* of the ordering. This
guards against wasting count files when the script logic is updated mid-run.

## Composition of the 5K seed set

```
1,383 unique sequence types (STs)
2,112 carbapenemase carriers (42.2%)

Klebsiella pneumoniae                            False  2,443
                                                 True   1,903
Klebsiella variicola subsp. variicola            False    231
                                                 True      72
Klebsiella quasipneumoniae subsp. similipneum.   False    135
                                                 True      94
Klebsiella quasipneumoniae subsp. quasipneum.    False     75
                                                 True      43
Klebsiella africana                              False      3
Klebsiella variicola subsp. tropica              False      1
```

The full per-tier breakdown is logged by the selection script and committed
to `assets/kpsc_expansion_subset_5k_meta.tsv` for the seed set.
