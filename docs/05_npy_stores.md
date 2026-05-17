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
| Plasmid family | 7 (one per AMR gene family) | sum 1 kb bins overlapping all member-variant gene regions in each family |

**Why gene families, not individual genes?** Nearly-identical allele variants
(NDM-1 vs NDM-5, OXA-48 vs OXA-181, CTX-M-15/14/65/27) sit on different
reference contigs. A read from a sample carrying any one variant maps equally
well to all of them; BWA places it at one position (MQ=0). Aggregating
across all family members sums the bin counts back into a single per-family
PCN signal — which also happens to match clinical reporting conventions
(AMRFinder+, Kleborate report at the family level by default).

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

## Plasmid-family NPY store

`data/setup/plasmid_genes_to_npy_kpsc.py` reads `plasmid_gene_coords.tsv` and
**sums** the 1 kb bins overlapping each gene region, then (with `--families`)
aggregates by gene family per `plasmid_gene_families.tsv`:

| Family | Member variants |
|---|---|
| `blaKPC` | blaKPC-2 |
| `blaCTX-M` | blaCTX-M-15, blaCTX-M-14, blaCTX-M-65, blaCTX-M-27 |
| `blaNDM` | blaNDM-1, blaNDM-5 |
| `blaOXA-48-like` | blaOXA-48, blaOXA-181 |
| `blaTEM` | blaTEM-1 |
| `qnrB` | qnrB1 |
| `aac6-Ib-cr` | aac6-Ib-cr |

The sum formula:

```python
count[sample, gene] = Σ bin.count
                      where bin.CHROM == gene.contig
                        AND bin.END    >  gene.start
                        AND bin.START  <  gene.end
```

CLI:

```bash
python3 data/setup/plasmid_genes_to_npy_kpsc.py \
    --counts-dir   data/raw/readcounts_subset_mq0 \
    --manifest     assets/kpsc_expansion_subset_5k.tsv \
    --gene-coords  assets/plasmid_refs/plasmid_gene_coords.tsv \
    --families     assets/plasmid_refs/plasmid_gene_families.tsv \
    --out-dir      data/inputs/KpSC-expansion-5k-mq0-plasmid-1000bp-npy
```

Output is `(n_samples, 7)` — one column per gene family. Downstream callers
and evaluation code use the family names (`blaKPC`, `blaCTX-M`, `blaNDM`,
`blaOXA-48-like`, `blaTEM`, `qnrB`, `aac6-Ib-cr`) directly.

## Build all four tiers

```bash
for N in 5k 10k 20k 40k; do
    # Chromosome store
    python3 data/setup/readcounts_to_npy_kpsc.py \
        --counts-dir data/raw/readcounts_subset_mq0 \
        --manifest   assets/kpsc_expansion_subset_${N}.tsv \
        --out-dir    data/inputs/KpSC-expansion-${N}-mq0-1000bp-npy \
        --keep-contigs NC_016845.1

    # Plasmid-family store
    python3 data/setup/plasmid_genes_to_npy_kpsc.py \
        --counts-dir  data/raw/readcounts_subset_mq0 \
        --manifest    assets/kpsc_expansion_subset_${N}.tsv \
        --gene-coords assets/plasmid_refs/plasmid_gene_coords.tsv \
        --families    assets/plasmid_refs/plasmid_gene_families.tsv \
        --out-dir     data/inputs/KpSC-expansion-${N}-mq0-plasmid-1000bp-npy
done
```

Each builder takes **5–15 minutes** depending on N and number of workers (the
TSV → numpy conversion is parallelised across `--workers` threads, default
16).
