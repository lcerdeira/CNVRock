"""
Download AllTheBacteria KpSC assemblies for Kleborate ground truth generation.

Strategy
--------
1. Download two small ATB metadata files (cached after first run):
     ena_metadata.tsv.gz        run_accession → sample_accession
     sample2species2file.tsv.gz sample_accession → tar.xz filename
2. For each run accession in kpsc_sra_accessions.txt, resolve:
     run_accession → sample_accession → klebsiella_*.asm.tar.xz
3. Download only the tar.xz files that contain at least one of our samples.
4. Extract only the FASTA files for our samples; rename to <run_accession>.fasta.

KpSC tar.xz files (ATB v0.2):
  klebsiella_pneumoniae__01-15.asm.tar.xz   (~15 files, large)
  klebsiella_quasipneumoniae__01.asm.tar.xz
  klebsiella_variicola__01.asm.tar.xz

Usage:
    # Limit to samples with completed BAMs (recommended):
    ls data/raw/bam/*.bam | xargs -I{} basename {} .bam > assets/kpsc_bam_accessions.txt

    python data/setup/download_assemblies.py \\
        --accessions   assets/kpsc_bam_accessions.txt \\
        --out-dir      data/assemblies/ \\
        --tmp-dir      data/assemblies/tmp/

    # Dry-run — print which tar.xz files would be downloaded:
    python data/setup/download_assemblies.py \\
        --accessions   assets/kpsc_bam_accessions.txt \\
        --dry-run
"""

import argparse
import gzip
import lzma
import os
import sys
import tarfile
import tempfile
import urllib.request


ATB_BASE = "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2"
ENA_META_URL    = f"{ATB_BASE}/metadata/ena_metadata.tsv.gz"
S2S_URL         = f"{ATB_BASE}/metadata/sample2species2file.tsv.gz"
ASSEMBLY_BASE   = f"{ATB_BASE}/assembly"

ENA_META_CACHE  = "assets/atb_ena_metadata.tsv.gz"
S2S_CACHE       = "assets/atb_sample2species2file.tsv.gz"

KPSC_SPECIES = {"klebsiella_pneumoniae", "klebsiella_quasipneumoniae",
                "klebsiella_variicola", "klebsiella_grimontii",
                "klebsiella_michiganensis", "klebsiella_oxytoca"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _download(url: str, dest: str, label: str = "") -> None:
    print(f"  Downloading {label or url} …", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    urllib.request.urlretrieve(url, dest)
    size_mb = os.path.getsize(dest) / 1024**2
    print(f"  → {dest}  ({size_mb:.1f} MB)", flush=True)


def load_ena_metadata(cache: str) -> dict[str, str]:
    """Return {run_accession: sample_accession}."""
    if not os.path.exists(cache):
        _download(ENA_META_URL, cache, "ena_metadata.tsv.gz")
    print(f"Parsing {cache} …", flush=True)
    run_to_sample: dict[str, str] = {}
    with gzip.open(cache, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            run_col    = header.index("run_accession")
            sample_col = header.index("sample_accession")
        except ValueError:
            raise RuntimeError(f"Expected columns not found. Header: {header[:10]}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(run_col, sample_col):
                run  = parts[run_col].strip()
                samp = parts[sample_col].strip()
                if run and samp:
                    run_to_sample[run] = samp
    print(f"  {len(run_to_sample):,} run→sample mappings loaded.", flush=True)
    return run_to_sample


def load_sample2species(cache: str) -> dict[str, str]:
    """Return {sample_accession: tar_filename}."""
    if not os.path.exists(cache):
        _download(S2S_URL, cache, "sample2species2file.tsv.gz")
    print(f"Parsing {cache} …", flush=True)
    sample_to_tar: dict[str, str] = {}
    with gzip.open(cache, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        try:
            samp_col = header.index("sample")
            spec_col = header.index("species")
            file_col = header.index("file")
        except ValueError:
            raise RuntimeError(f"Expected columns not found. Header: {header}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(samp_col, spec_col, file_col):
                samp = parts[samp_col].strip()
                spec = parts[spec_col].strip()
                tar  = parts[file_col].strip()
                if samp and spec in KPSC_SPECIES and tar:
                    sample_to_tar[samp] = tar
    print(f"  {len(sample_to_tar):,} KpSC sample→tar mappings loaded.", flush=True)
    return sample_to_tar


# ---------------------------------------------------------------------------
# Extract FASTAs from a tar.xz
# ---------------------------------------------------------------------------

def _extract_from_tarxz(tar_path: str, wanted_samples: set[str],
                         run_for_sample: dict[str, str], out_dir: str) -> dict[str, str]:
    """Extract FASTAs for wanted_samples from tar_path.

    Returns {run_accession: fasta_path} for successfully extracted samples.
    """
    extracted: dict[str, str] = {}
    print(f"  Extracting from {os.path.basename(tar_path)} …", flush=True)

    with tarfile.open(tar_path, "r:xz") as tf:
        members = tf.getmembers()
        for member in members:
            name = os.path.basename(member.name)
            # FASTAs inside tar are named like SAMEA12345.fasta or SAMEA12345.fa
            stem = name.split(".")[0]
            if stem not in wanted_samples:
                continue
            run_acc = run_for_sample[stem]
            dest    = os.path.join(out_dir, f"{run_acc}.fasta")
            if os.path.exists(dest) and os.path.getsize(dest) > 500:
                extracted[run_acc] = dest
                continue
            f = tf.extractfile(member)
            if f is None:
                continue
            data = f.read()
            # Decompress if gzipped inside the tar
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            with open(dest, "wb") as out:
                out.write(data)
            extracted[run_acc] = dest
            print(f"    extracted {run_acc} ({len(data)//1024} KB)", flush=True)

    return extracted


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accessions",  default="assets/kpsc_bam_accessions.txt",
                        help="One SRA run accession per line (default: kpsc_bam_accessions.txt).")
    parser.add_argument("--out-dir",     default="data/assemblies/",
                        help="Output directory for extracted FASTA files.")
    parser.add_argument("--tmp-dir",     default=None,
                        help="Temp directory for tar.xz downloads (default: system tmp).")
    parser.add_argument("--dry-run",     action="store_true",
                        help="Print plan without downloading anything.")
    args = parser.parse_args()

    with open(args.accessions) as f:
        run_accessions = [l.strip() for l in f if l.strip()]
    print(f"Run accessions to resolve: {len(run_accessions)}", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Build run → sample → tar mapping ─────────────────────────────────
    run_to_sample  = load_ena_metadata(ENA_META_CACHE)
    sample_to_tar  = load_sample2species(S2S_CACHE)

    # For each run accession, resolve sample + tar
    tar_to_runs: dict[str, list[tuple[str, str]]] = {}  # tar → [(run, sample)]
    unresolved: list[str] = []

    for run in run_accessions:
        sample = run_to_sample.get(run)
        if sample is None:
            unresolved.append(run)
            continue
        tar = sample_to_tar.get(sample)
        if tar is None:
            unresolved.append(run)
            continue
        tar_to_runs.setdefault(tar, []).append((run, sample))

    print(f"\nResolved:   {sum(len(v) for v in tar_to_runs.values())} samples across "
          f"{len(tar_to_runs)} tar.xz files", flush=True)
    print(f"Unresolved: {len(unresolved)} run accessions (not in ATB KpSC)", flush=True)

    if unresolved:
        print(f"  First 10 unresolved: {unresolved[:10]}", flush=True)

    print("\nTar files needed:", flush=True)
    for tar, runs in sorted(tar_to_runs.items()):
        # Count already done
        n_done = sum(
            1 for run, _ in runs
            if os.path.exists(os.path.join(args.out_dir, f"{run}.fasta"))
            and os.path.getsize(os.path.join(args.out_dir, f"{run}.fasta")) > 500
        )
        print(f"  {tar:<50}  {len(runs)} samples  ({n_done} already done)", flush=True)

    if args.dry_run:
        return

    # ── Download and extract each required tar.xz ─────────────────────────
    tmp_root = args.tmp_dir or tempfile.gettempdir()
    os.makedirs(tmp_root, exist_ok=True)

    total_ok, total_skip, total_fail = 0, 0, 0

    for tar_name, runs in sorted(tar_to_runs.items()):
        # Skip if all samples already extracted
        pending = [
            (run, samp) for run, samp in runs
            if not (os.path.exists(os.path.join(args.out_dir, f"{run}.fasta"))
                    and os.path.getsize(os.path.join(args.out_dir, f"{run}.fasta")) > 500)
        ]
        already_done = len(runs) - len(pending)
        total_skip += already_done
        if not pending:
            print(f"\n{tar_name}: all {len(runs)} already extracted — skipping.", flush=True)
            continue

        url      = f"{ASSEMBLY_BASE}/{tar_name}"
        tar_path = os.path.join(tmp_root, tar_name)

        print(f"\n{'='*60}", flush=True)
        print(f"Processing {tar_name} ({len(pending)} samples needed) …", flush=True)

        if not os.path.exists(tar_path):
            try:
                _download(url, tar_path, tar_name)
            except Exception as e:
                print(f"  ERROR downloading {tar_name}: {e}", file=sys.stderr)
                total_fail += len(pending)
                continue

        wanted   = {samp for _, samp in pending}
        run_map  = {samp: run for run, samp in pending}
        try:
            extracted = _extract_from_tarxz(tar_path, wanted, run_map, args.out_dir)
            total_ok += len(extracted)
            missing = wanted - set(run_map[r] for r in extracted)  # type: ignore[arg-type]
            if missing:
                print(f"  WARNING: {len(missing)} samples not found in tar: "
                      f"{list(missing)[:5]}", file=sys.stderr)
                total_fail += len(missing)
        except Exception as e:
            print(f"  ERROR extracting {tar_name}: {e}", file=sys.stderr)
            total_fail += len(pending)

        # Remove tar after extraction to save disk space
        try:
            os.unlink(tar_path)
            print(f"  Removed {tar_path}", flush=True)
        except OSError:
            pass

    print(f"\n{'='*60}", flush=True)
    print(f"Done.  extracted={total_ok}  skipped={total_skip}  failed={total_fail}",
          flush=True)
    n_fasta = len([f for f in os.listdir(args.out_dir) if f.endswith(".fasta")])
    print(f"Assemblies in {args.out_dir}: {n_fasta}", flush=True)


if __name__ == "__main__":
    main()
