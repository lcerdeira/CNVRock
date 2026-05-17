# 2. Data acquisition

KpSC short-read sequencing data lives on ENA, addressed by run accession
(`ERR…`, `DRR…`, `SRR…`). For ~78,000 samples this is **multiple TB of FASTQ**,
which precludes brittle HTTPS streaming.

## Why Aspera, not wget

We tried both:

| Mode | Single-task | 10 concurrent | 50 concurrent |
|---|---|---|---|
| `wget` (HTTPS) | ✅ works | ❌ R2 fails — instant rejection | ❌ all fail |
| `ascp` (Aspera) | ✅ works | ✅ 10/10 succeed in ~46s | ❌ ~50% fail |

EBI imposes an undocumented per-source-IP **concurrent connection limit** on
HTTPS that drops R2 (the second connection per task) almost instantly. Aspera's
FASP protocol multiplexes over a single SSH session and is the
**EBI-recommended bulk-download protocol**, so we use it throughout.

For details, see EBI's [Bio Image Archive download
guide](https://www.ebi.ac.uk/bioimage-archive/help-download/), which is where
the public-domain ed25519 SSH key used by `era-fasp@fasp.sra.ebi.ac.uk` is
distributed.

## Setup on HPC

```bash
# In the cnvrock conda env:
conda install -y -c bioconda aspera-cli

# The bioconda package ships ascli (Ruby wrapper) but not ascp itself.
# Install ascp via:
ascli config ascp install
# Drops the binary at ~/.aspera/sdk/ascp

# Save the EBI public key (one-time):
curl -sL https://www.ebi.ac.uk/bioimage-archive/help-download/ \
    | sed -n '/BEGIN OPENSSH PRIVATE KEY/,/END OPENSSH PRIVATE KEY/p' \
    > ~/.aspera/sdk/ebi_aspera_key.openssh
chmod 600 ~/.aspera/sdk/ebi_aspera_key.openssh
```

## Per-task pipeline

`hpc/aspera_subset_pipeline.sh` runs per SLURM array task:

1. Read one row of the **subset manifest** TSV (cols: `accession`, `layout`,
   `r1_url`, `r2_url`)
2. **`ascp` download** R1 (and R2 if PAIRED) to `data/raw/fastq_subset/`
3. **BWA-MEM** paired-end align to `HS11286_extended.fasta`, then `samtools
   sort` and `samtools index`. Output to `data/raw/bam_subset/`.
4. **GATK CollectReadCounts** at 1 kb resolution across chrom + plasmid
   contigs, with `--minimum-mapping-quality 20`. Output count TSV to
   `data/raw/readcounts_subset_mq20/{ACC}.counts.tsv`.

FASTQs and BAMs are **kept on disk** so we can re-extract counts at different
MQ thresholds without re-downloading.

## SLURM submission

Concurrency is capped at **10 per array** to respect Aspera's per-IP limit:

```bash
# Submit the 20K manifest as four chunks of ~5,000 array tasks each
# (SLURM MaxArraySize=5000 on LSHTM HPC).
for chunk in 1 2 3 4; do
    case $chunk in
        1) OFFSET=0;     ARR=1-4999 ;;
        2) OFFSET=4999;  ARR=1-4999 ;;
        3) OFFSET=9998;  ARR=1-4999 ;;
        4) OFFSET=14997; ARR=1-5000 ;;
    esac
    sbatch --parsable \
        --export=ALL,MANIFEST=$REPO/assets/kpsc_expansion_subset_20k.tsv,BATCH_OFFSET=$OFFSET \
        --array=$ARR%10 hpc/aspera_subset_pipeline.sh
done
```

## Throughput

| Stage | Wall time per sample (avg) |
|---|---|
| `ascp` download (R1 + R2, ~300 MB total) | 30–60 s |
| BWA mem align (4 CPUs, ~5 M reads, 7.2 Mb ref) | 2–3 min |
| `samtools sort` + `index` | ~30 s |
| GATK CollectReadCounts | ~30 s |
| **Total per task** | **~3–5 min median** |

At concurrency 10, the pipeline sustains ~150–250 samples/hour. The 5K tier
took ~30 hours; the 10K and 20K tiers add ~25 and ~50 hours respectively.

## Compliance with HPC network policy

A single-IP HTTPS surge from the HPC earlier caused a brief outage. We
coordinated with HPC support to limit concurrent connections; Aspera was
explicitly endorsed as a low-risk alternative. The CONCUR=10 setting is well
within EBI's published guidance.
