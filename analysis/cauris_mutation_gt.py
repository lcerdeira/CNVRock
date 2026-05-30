#!/usr/bin/env python3
"""
Build the C. auris ERG11 azole-resistance mutation ground truth.

For the 500-isolate subset re-aligned to B8441v3 (BAMs in bam_cauris_sub500),
call the canonical ERG11 azole-resistance substitutions (Y132F, K143R, F126L)
per isolate via a bcftools consensus → translate → compare-to-reference
approach. The resulting per-sample R-mutation flag lets us test the
manuscript's secondary claim: does CNVRock-detected ERG11 *copy-number
amplification* predict azole resistance independently of these point
mutations?

ERG11 / CYP51 (lanosterol 14-alpha demethylase): B8441 CM076440.1,
1,474,722-1,476,296, minus strand (locus B9J08_03698).

Run on HPC (samtools + bcftools in cnvrock env), after the 500 BAMs exist:
  python3 analysis/cauris_mutation_gt.py
"""
from __future__ import annotations
import os, glob, subprocess, tempfile
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

REPO = Path("/home/lshlt19/CNVRock")
REF = REPO / "assets/cauris_ref/B8441v3.fasta"
BAM_DIR = REPO / "data/raw/bam_cauris_sub500"
OUT = REPO / "data/results/cauris_mutation_gt"
GENE_CALLS = REPO / "data/results/39_cauris/gene_calls.tsv"
CONTIG, START, END, STRAND = "CM076440.1", 1474722, 1476296, "-"
REGION = f"{CONTIG}:{START}-{END}"
# canonical ERG11 azole-resistance substitutions in C. auris
HOTSPOTS = {126: ("F", "L"), 132: ("Y", "F"), 143: ("K", "R")}

# samtools path — the cnvrock env (Biopython) doesn't ship samtools;
# use the aligners env which has it. Fall back to plain "samtools" if
# the explicit path doesn't exist (e.g. on a different cluster node).
_SAMTOOLS_EXPLICIT = "/home/lshlt19/miniconda3/envs/aligners/bin/samtools"
SAMTOOLS = _SAMTOOLS_EXPLICIT if Path(_SAMTOOLS_EXPLICIT).exists() else "samtools"


def ref_protein():
    for rec in SeqIO.parse(str(REF), "fasta"):
        if rec.id == CONTIG:
            cds = rec.seq[START - 1:END]
            if STRAND == "-":
                cds = cds.reverse_complement()
            return cds.translate()
    raise RuntimeError("contig not found")


def sample_protein(bam, scratch):
    """Consensus CDS for one BAM over the ERG11 region, translated."""
    cons = os.path.join(scratch, "cons.fa")
    # samtools consensus restricted to the gene region
    with open(cons, "w") as fh:
        subprocess.run([SAMTOOLS, "consensus", "-r", REGION, "-f", "fasta", bam],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)
    recs = list(SeqIO.parse(cons, "fasta"))
    if not recs:
        return None
    seq = recs[0].seq
    # consensus is the full region on + strand; take CDS and orient
    if STRAND == "-":
        seq = seq.reverse_complement()
    # trim to multiple of 3
    seq = seq[: len(seq) - (len(seq) % 3)]
    try:
        return seq.translate()
    except Exception:
        return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    refp = ref_protein()
    print(f"ERG11 reference protein length {len(refp)}; "
          f"hotspot ref AAs: " + ", ".join(
              f"{p}{refp[p-1]}" for p in sorted(HOTSPOTS)))
    bams = sorted(glob.glob(str(BAM_DIR / "*.bam")))
    print(f"BAMs found: {len(bams)}")
    rows = []
    for bam in bams:
        acc = os.path.basename(bam)[:-4]
        with tempfile.TemporaryDirectory() as scratch:
            try:
                prot = sample_protein(bam, scratch)
            except Exception:
                prot = None
        muts = []
        if prot is not None and len(prot) >= max(HOTSPOTS):
            for pos, (wt, mut) in HOTSPOTS.items():
                aa = str(prot[pos - 1]) if pos - 1 < len(prot) else "?"
                if aa == mut:
                    muts.append(f"{wt}{pos}{mut}")
        rows.append({"accession": acc,
                     "erg11_R_mutation": int(len(muts) > 0),
                     "mutations": ";".join(muts) if muts else "-"})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "cauris_erg11_mutation_gt.tsv", sep="\t", index=False)
    n_mut = int(df["erg11_R_mutation"].sum())
    print(f"\nwrote {OUT/'cauris_erg11_mutation_gt.tsv'} ({len(df)} isolates)")
    print(f"  isolates with an ERG11 azole-R mutation: {n_mut}")

    # cross with CNVRock ERG11 copy number — the headline secondary analysis
    if GENE_CALLS.exists() and n_mut > 0:
        g = pd.read_csv(GENE_CALLS, sep="\t")[["sample_id", "crr_ERG11", "ERG11"]]
        m = df.merge(g, left_on="accession", right_on="sample_id", how="inner")
        amp = m["crr_ERG11"] > 1.5
        mut = m["erg11_R_mutation"] == 1
        print(f"\n  ERG11 status among {len(m)} isolates with both:")
        print(f"    mutation only (no amp):   {int((mut & ~amp).sum())}")
        print(f"    amplification only (WT):  {int((~mut & amp).sum())}  "
              f"<- CNV-only azole-R candidates")
        print(f"    both:                     {int((mut & amp).sum())}")
        print(f"    neither:                  {int((~mut & ~amp).sum())}")


if __name__ == "__main__":
    main()
