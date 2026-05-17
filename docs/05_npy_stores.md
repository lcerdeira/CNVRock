# 5. NPY stores

CNVRock's training pipeline reads from **numpy NPY stores**, not the raw GATK
count TSVs. Each store is a directory with:

```
data/inputs/KpSC-expansion-{N}-mq20-1000bp-npy/
├── contigs.npy      structured array [(chrom, start, end)] per bin
├── counts.npy       uint32 (n_samples, n_bins)
└── sample_ids.npy   object array of run accessions, matches counts row order
```

Per scaling tier we materialise **two stores**:

| Store | Bins | Source |
|---|---|---|
| Chromosome | 5,334 (NC_016845.1 only) | filter 1 kb bins by contig |
| Plasmid-gene | 12 (one per AMR gene) | sum 1 kb bins overlapping each gene region |

The combined `aspera_subset_pipeline.sh` produces a single count TSV per
sample with all 7,365 bins. Two separate builders then materialise the two
stores from those same files.

## Chromosome NPY store

`data/setup/readcounts_to_npy_kpsc.py` keeps every 1 kb bin in
`--keep-contigs`. For the chromosome store:

```bash
python3 data/setup/readcounts_to_npy_kpsc.py \
    --counts-dir   data/raw/readcounts_subset_mq20 \
    --manifest     assets/kpsc_expansion_subset_5k.tsv \
    --out-dir      data/inputs/KpSC-expansion-5k-mq20-1000bp-npy \
    --keep-contigs NC_016845.1
```

Validates that every sample has the same number of bins; raises if not.

## Plasmid-gene NPY store

`data/setup/plasmid_genes_to_npy_kpsc.py` reads `plasmid_gene_coords.tsv` and
**sums** the 1 kb bins overlapping each gene region:

```python
count[sample, gene] = Σ bin.count
                      where bin.CHROM == gene.contig
                        AND bin.END    >  gene.start
                        AND bin.START  <  gene.end
```

CLI:

```bash
python3 data/setup/plasmid_genes_to_npy_kpsc.py \
    --counts-dir   data/raw/readcounts_subset_mq20 \
    --manifest     assets/kpsc_expansion_subset_5k.tsv \
    --gene-coords  assets/plasmid_refs/plasmid_gene_coords.tsv \
    --out-dir      data/inputs/KpSC-expansion-5k-mq20-plasmid-1000bp-npy
```

Output is `(n_samples, 12)` — one column per gene of interest. Matches the
Phase 1 plasmid store schema, so downstream callers and evaluation code
require no changes.

## Build all four tiers

```bash
for N in 5k 10k 20k 40k; do
    # Chromosome store
    python3 data/setup/readcounts_to_npy_kpsc.py \
        --counts-dir data/raw/readcounts_subset_mq20 \
        --manifest   assets/kpsc_expansion_subset_${N}.tsv \
        --out-dir    data/inputs/KpSC-expansion-${N}-mq20-1000bp-npy \
        --keep-contigs NC_016845.1

    # Plasmid-gene store
    python3 data/setup/plasmid_genes_to_npy_kpsc.py \
        --counts-dir data/raw/readcounts_subset_mq20 \
        --manifest   assets/kpsc_expansion_subset_${N}.tsv \
        --gene-coords assets/plasmid_refs/plasmid_gene_coords.tsv \
        --out-dir    data/inputs/KpSC-expansion-${N}-mq20-plasmid-1000bp-npy
done
```

Each builder takes **5–15 minutes** depending on N and number of workers (the
TSV → numpy conversion is parallelised across `--workers` threads, default
16).
