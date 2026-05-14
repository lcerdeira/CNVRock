# Experiment 21 — KpSC Chromosomal CNV (first microbial application)

**Status:** Data preparation in progress

## Biological Question

Does the deep-learning VAE + HMM approach for CNV discovery generalise from
*Plasmodium falciparum* (malaria) to bacterial pathogens? Specifically: can
the pipeline detect chromosomal gene copy number variation in the
**Klebsiella pneumoniae Species Complex (KpSC)** — a WHO critical-priority
pathogen causing carbapenem-resistant nosocomial infections?

## Why KpSC

- Clinically urgent: carbapenem-resistant KpSC kills tens of thousands annually
- Well-studied genome: HS11286 (NC_016845.1) is a high-quality MDR reference
- Large public dataset: AllTheBacteria provides 100k+ KpSC assemblies with
  associated SRA accessions for raw reads
- Tractable biology: CNVs of chromosomal resistance genes are well-documented
  (blaSHV amplification, porin deletion)

## What's New vs Experiment 20 (Pf)

| Component | Exp 20 (Pf) | Exp 21 (KpSC) |
|---|---|---|
| Organism | *P. falciparum* | *K. pneumoniae* Species Complex |
| Reference genome | PlasmoDB-54 Pf3D7 (20.8 Mb, 16 contigs) | HS11286 NC_016845.1 (~5.3 Mb, 1 chromosome) |
| Bins | 20,814 (core genome) | ~4,000–5,000 (core genome, TBD after Panaroo) |
| VAE architecture | 05_conv_vae (Pf constants) | 06_conv_vae (n_bins_raw from data) |
| Genes of interest | MDR1, CRT, GCH1, PM2_PM3 | blaSHV, ompK35, ompK36, ramR |
| Ground truth | Pf8 GATK CNV calls | AMRFinder+ on AllTheBacteria assemblies |
| Data source | MalariaGEN Pf9 (53,973 samples) | NCBI SRA via AllTheBacteria catalog |
| CNV downsampling | Yes (pf9_gt_path) | No (not yet — GT being built) |

## Data Preparation Steps (required before training)

1. **Sample selection** (`data/setup/atb_sample_selection.py`):
   - Download AllTheBacteria KpSC metadata
   - Filter: Illumina, QC-pass, estimated coverage ≥ 30×
   - Extract SRA accessions → `assets/kpsc_sra_accessions.txt`

2. **Raw read download** (NCBI SRA toolkit):
   ```bash
   while read acc; do
     fasterq-dump "$acc" --outdir fastq/ --split-files
   done < assets/kpsc_sra_accessions.txt
   ```

3. **Map to HS11286 reference** (bwa mem):
   ```bash
   # Download reference first:
   # ncbi-datasets download genome accession GCF_000240185.1 --include genome,gff3
   bwa mem -t 8 HS11286.fasta sample_R1.fastq.gz sample_R2.fastq.gz \
     | samtools sort -o sample.bam && samtools index sample.bam
   ```

4. **Extract gene coordinates** from HS11286 GFF (needed for 06_gene_cnv_caller.py):
   ```bash
   grep -E "blaSHV|ompK35|ompK36|ramR" HS11286.gff3 | \
     awk '{print $1, $4, $5, $9}' > gene_coords.txt
   # Then update GENES_OF_INTEREST in models/cnv/06_gene_cnv_caller.py
   ```

5. **Core genome BED** (using Panaroo on AllTheBacteria assemblies):
   ```bash
   panaroo -i assemblies/*.gff -o panaroo_out/ --clean-mode strict
   # Convert core gene intervals → core-genome.bed for HS11286 coordinate space
   ```

6. **Run GATK CollectReadCounts** (Nextflow — see main.nf in this directory):
   ```bash
   nextflow run main.nf -profile local -resume
   ```

7. **Convert TSV → NPY** (readcounts_to_npy.py in this directory):
   ```bash
   .venv/bin/python readcounts_to_npy.py
   ```

8. **Generate AMRFinder+ ground truth**:
   ```bash
   # Run on AllTheBacteria assemblies for samples with SRA reads
   for asm in assemblies/*.fasta; do
     sample=$(basename "$asm" .fasta)
     amrfinder --plus -O Klebsiella_pneumoniae -n "$asm" --output gt/${sample}.tsv
   done
   # Then run data/setup/atb_sample_selection.py --summarise-amrfinder to aggregate
   ```

9. **Train** (once data is ready):
   ```bash
   cd /path/to/deep-gw-cnv/models
   .venv/bin/python train.py experiments/21/config.yaml
   ```

## Genes of Interest

| Gene | Type | AMRFinder+ event | Clinical significance |
|---|---|---|---|
| blaSHV | amp | count ≥ 2 | Extra chromosomal copies → elevated beta-lactam MIC |
| ompK35 | del | count = 0 | Porin loss → reduced drug uptake → carbapenem resistance |
| ompK36 | del | count = 0 | Porin loss (OmpC homolog) → same |
| ramR   | del | count = 0 | Efflux repressor loss → AcrAB-TolC upregulation → MDR |

## Hypothesis

The VAE will learn that the typical KpSC sample has CN=1 across the core
chromosome. Amplification of blaSHV (CN≥2) and deletion of porin/ramR genes
(CN=0) will appear as anomalies in the latent space, producing elevated
reconstruction error at those loci, which the HMM will segment as CN≠1 states.

Expected challenges vs Pf:
- Higher within-species genetic diversity → VAE reconstruction harder
- Fewer bins (5k vs 21k) → smaller model, potentially faster convergence
- No equivalent of "AF-E historical" population bias
- Ground truth is assembly-derived (not validated by independent pipeline)
