# 3. Reference and intervals

## `HS11286_extended.fasta`

The base reference is the **HS11286** assembly of *K. pneumoniae*, augmented
with **11 plasmid contigs** covering the Phase D AMR panel:

| Contig | Length | Genes of interest |
|---|---|---|
| NC_016845.1 | 5,333,942 bp | KpSC chromosome (carries `blaSHV`) |
| NC_016846.1 | 111,195 bp | `blaKPC-2` |
| MK552109.1 | 180,286 bp | `blaCTX-M-15` |
| MZ606384.2 | 140,133 bp | `blaNDM-1` |
| NZ_CP031850.1 | 239,020 bp | `blaCTX-M-14` |
| NZ_KX236178.1 | 121,348 bp | `blaTEM-1` |
| NC_016980.1 | 267,242 bp | `qnrB1` |
| CP006662.2 | 85,161 bp | `aac(6')-Ib-cr` |
| JN626286.1 | 61,881 bp | `blaOXA-48` |
| CP034201.2 | 372,826 bp | `blaNDM-5` |
| CP113224.1 | 53,691 bp | `blaOXA-181` |
| CP030319.1 | 169,824 bp | `blaCTX-M-65` |
| CP138680.1 | 221,802 bp | `blaCTX-M-27` |

**Total: ~7.2 Mb** (chromosome + 11 plasmid backbones).

## 1 kb interval list

CNVRock trains on per-bin read counts. The interval list defines those bins;
**this used to be a bug source** — the original Phase 2 pipeline pointed at
`kpsc-whole-chrom.interval_list` which declared the chromosome as a *single*
interval. GATK `CollectReadCounts` does **not subdivide** intervals, so it
produced one-bin-per-sample count files unusable for CNVRock.

Fixed by regenerating with explicit 1 kb binning:

```bash
gatk PreprocessIntervals \
    --reference HS11286_extended.fasta \
    --bin-length 1000 \
    --padding 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output HS11286_extended_1kb.interval_list
```

Output: **7,365 bins** total = 5,334 on NC_016845.1 + 2,031 across the 12
other contigs.

## Plasmid gene coordinates

`assets/plasmid_refs/plasmid_gene_coords.tsv` defines the 12 plasmid gene
regions used for **gene-level PCN aggregation**:

```
gene            contig         start    end     absent_threshold
blaKPC-2        NC_016846.1    20557    21435   0.20
blaCTX-M-15     MK552109.1     119392   120264  0.50
blaNDM-1        MZ606384.2     90937    91746   0.20
blaCTX-M-14     NZ_CP031850.1  113509   114381  1.00
blaTEM-1        NZ_KX236178.1  66790    67650   0.20
qnrB1           NC_016980.1    18421    19100   0.20
aac6-Ib-cr      CP006662.2     38208    38689   0.10
blaOXA-48       JN626286.1     5445     6242    0.20
blaNDM-5        CP034201.2     69698    70263   0.20
blaOXA-181      CP113224.1     28840    29637   0.20
blaCTX-M-65     CP030319.1     6227     7049    0.20
blaCTX-M-27     CP138680.1     183201   184385  0.20
```

The `absent_threshold` is the **per-gene PCN cutoff** below which the plasmid
is considered absent (handles low-mappability/duplicate-sequence cases).

## Mapping quality threshold

GATK's `--minimum-mapping-quality 20` keeps reads where BWA estimates ≤1%
mis-placement probability. See {doc}`09_methods` for the rationale (the
standard GATK HaplotypeCaller / DELLY default), and why we **lowered from
MQ=40** (the Phase 1 value) — the extended-reference plasmids share sequence
with each other, so MQ=40 was discarding the multi-mapping reads we need for
plasmid gene-level CNV calls.
