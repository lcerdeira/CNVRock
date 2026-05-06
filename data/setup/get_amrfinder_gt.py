"""
Download KpSC assemblies from AllTheBacteria and run AMRFinder+ to generate
ground truth gene copy counts for evaluation.

Streams each species chunk from ATB EBI FTP, extracts only the assemblies
for our target samples (no full chunk saved to disk), runs AMRFinder+
immediately, then deletes the FASTA. Writes:

    assets/kpsc_amrfinder_ground_truth.tsv
        Columns: sample_id, blaSHV, ompK35, ompK36, ramR
        (sample_id = run accession, copy count = integer)

Usage (run from CNVRock/ repo root):
    python3 data/setup/get_amrfinder_gt.py
"""

import csv
import os
import subprocess
import tarfile
import tempfile
import time
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

REPO_DIR      = Path("/home/lshlt19/CNVRock")
ASSETS        = REPO_DIR / "assets"
TMP_DIR       = REPO_DIR / "data" / "amrfinder_tmp"
OUT_TSV       = ASSETS   / "kpsc_amrfinder_ground_truth.tsv"
RAW_DIR       = REPO_DIR / "data" / "amrfinder_raw"   # one .tsv per sample

ATB_BASE      = "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly"

# All Klebsiella species chunks to search (covers KpSC: Kp, Kv, Kq, Ka, etc.)
SPECIES_PREFIXES = [
    "klebsiella_pneumoniae",
    "klebsiella_variicola",
    "klebsiella_quasipneumoniae",
    "klebsiella_aerogenes",
    "klebsiella_michiganensis",
    "klebsiella_oxytoca",
    "klebsiella_grimontii",
    "klebsiella_ornithinolytica",
    "klebsiella_planticola",
]

# AMRFinder gene names to count (partial match on Gene symbol column).
# Chromosomal genes: blaSHV, ompK35, ompK36, ramR
# Plasmid genes:     blaKPC, blaCTX-M, blaNDM  (any variant in the family counts)
GENE_PATTERNS = {
    "blaSHV":    ["blaSHV"],
    "ompK35":    ["ompK35", "OmpK35"],
    "ompK36":    ["ompK36", "OmpK36"],
    "ramR":      ["ramR",   "RamR"],
    "blaKPC":    ["blaKPC"],
    "blaCTX-M":  ["blaCTX-M"],
    "blaNDM":    ["blaNDM"],
    # Phase C plasmid genes (added for exp 28 evaluation)
    "qnrB1":     ["qnrB1"],
    "blaOXA-48": ["OXA-48", "blaOXA-48"],
    "aac6-Ib-cr":["aac(6')-Ib-cr", "aac(6)-Ib-cr"],
}

# ---------------------------------------------------------------------------
# Step 1: build run→sample and sample→run maps
# ---------------------------------------------------------------------------

def build_sample_map():
    run_to_sample = {}
    with open(ASSETS / "kpsc_sample_run_map.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            run_to_sample[row["run_accession"]] = row["sample_accession"]

    run_ids = []
    with open(ASSETS / "kpsc-paths-to-bams.tsv") as f:
        for line in f:
            rid = line.split("\t")[0].strip()
            if rid:
                run_ids.append(rid)

    sample_to_run = {}
    for rid in run_ids:
        if rid in run_to_sample:
            sa = run_to_sample[rid]
            sample_to_run[sa] = rid

    print(f"Target samples: {len(sample_to_run)} (from {len(run_ids)} run IDs)")
    return sample_to_run


# ---------------------------------------------------------------------------
# Step 2: stream each ATB chunk, extract matching assemblies, run AMRFinder
# ---------------------------------------------------------------------------

def run_amrfinder(fasta_path: Path, sample_id: str, out_path: Path) -> bool:
    """Run AMRFinder+ on a single FASTA. Returns True on success."""
    cmd = [
        "amrfinder",
        "-O", "Klebsiella_pneumoniae",
        "-n", str(fasta_path),
        "--output", str(out_path),
        "--quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  AMRFinder failed for {sample_id}: {result.stderr[:200]}")
        return False
    return True


def process_chunk(url: str, needed: set, sample_to_run: dict,
                  raw_dir: Path, already_done: set) -> set:
    """Stream one chunk, extract and process assemblies for `needed` samples.
    Returns set of sample accessions successfully processed."""
    remaining = needed - already_done
    if not remaining:
        return set()

    print(f"  Streaming {url.split('/')[-1]} ...", flush=True)
    t0 = time.time()

    proc = subprocess.Popen(["curl", "-sf", url], stdout=subprocess.PIPE)
    xz   = subprocess.Popen(["xz", "-d"], stdin=proc.stdout, stdout=subprocess.PIPE)

    done_this_chunk = set()

    try:
        with tarfile.open(fileobj=xz.stdout, mode="r|") as tf:
            for member in tf:
                name   = Path(member.name).name          # e.g. SAMN24219486.fa
                sample = name.replace(".fa", "").replace(".fasta", "")
                if sample not in remaining:
                    continue

                run_id   = sample_to_run[sample]
                out_path = raw_dir / f"{run_id}.amrfinder.tsv"
                if out_path.exists():
                    done_this_chunk.add(sample)
                    continue

                # Extract to temp file
                with tempfile.NamedTemporaryFile(suffix=".fa", delete=False,
                                                 dir=TMP_DIR) as tmp:
                    tmp_path = Path(tmp.name)
                    fobj = tf.extractfile(member)
                    if fobj is None:
                        continue
                    tmp.write(fobj.read())

                ok = run_amrfinder(tmp_path, run_id, out_path)
                tmp_path.unlink(missing_ok=True)

                if ok:
                    done_this_chunk.add(sample)
                    print(f"    {run_id} ({sample}) done  [{len(done_this_chunk)} this chunk]",
                          flush=True)

    except Exception as exc:
        print(f"  Error processing chunk: {exc}")
    finally:
        proc.kill()
        xz.kill()

    elapsed = time.time() - t0
    print(f"  Chunk done in {elapsed:.0f}s, found {len(done_this_chunk)} samples", flush=True)
    return done_this_chunk


# ---------------------------------------------------------------------------
# Step 3: parse AMRFinder TSVs → gene copy counts
# ---------------------------------------------------------------------------

def count_gene_copies(tsv_path: Path) -> dict:
    counts = {g: 0 for g in GENE_PATTERNS}
    if not tsv_path.exists():
        return counts
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_sym = row.get("Gene symbol", "") or row.get("Gene Symbol", "")
            for gene, patterns in GENE_PATTERNS.items():
                if any(p in gene_sym for p in patterns):
                    counts[gene] += 1
    return counts


def aggregate_results(sample_to_run: dict, raw_dir: Path) -> None:
    rows = []
    for sample, run_id in sample_to_run.items():
        tsv = raw_dir / f"{run_id}.amrfinder.tsv"
        counts = count_gene_copies(tsv)
        rows.append({"sample_id": run_id, **counts})

    rows.sort(key=lambda r: r["sample_id"])
    genes = list(GENE_PATTERNS.keys())  # blaSHV, ompK35, ompK36, ramR, blaKPC, blaCTX-M, blaNDM

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id"] + genes, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    present = sum(1 for r in rows if any(r[g] > 0 for g in genes))
    print(f"\nWrote {len(rows)} samples to {OUT_TSV}")
    print(f"  Samples with ≥1 gene detected: {present}/{len(rows)}")
    for g in genes:
        n_amp = sum(1 for r in rows if r[g] >= 2)
        n_one = sum(1 for r in rows if r[g] == 1)
        n_zero = sum(1 for r in rows if r[g] == 0)
        n_missing = sum(1 for r in rows if not (raw_dir / f"{r['sample_id']}.amrfinder.tsv").exists())
        print(f"  {g}: 0={n_zero}  1={n_one}  ≥2={n_amp}  no_result={n_missing}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    sample_to_run = build_sample_map()
    needed        = set(sample_to_run.keys())

    # Check what's already done (resume-safe)
    already_done = set()
    for sample, run_id in sample_to_run.items():
        if (RAW_DIR / f"{run_id}.amrfinder.tsv").exists():
            already_done.add(sample)
    if already_done:
        print(f"Resuming: {len(already_done)} already done, {len(needed)-len(already_done)} remaining")

    found = set(already_done)

    for sp in SPECIES_PREFIXES:
        for chunk_n in range(1, 50):
            url   = f"{ATB_BASE}/{sp}__{chunk_n:02d}.asm.tar.xz"
            check = subprocess.run(
                ["curl", "-sf", "--head", url],
                capture_output=True, text=True
            )
            if "404" in check.stdout or "404" in check.stderr or check.returncode != 0:
                break   # no more chunks for this species

            done = process_chunk(url, needed, sample_to_run, RAW_DIR, found)
            found.update(done)

            if found >= needed:
                print("All samples found!")
                break
        if found >= needed:
            break

    missing = needed - found
    if missing:
        runs = [sample_to_run[s] for s in missing]
        print(f"\n{len(missing)} samples not found in any ATB chunk:")
        for r in runs[:20]:
            print(f"  {r}")

    print("\nAggregating AMRFinder results...")
    aggregate_results(sample_to_run, RAW_DIR)


if __name__ == "__main__":
    main()
