"""
Download AllTheBacteria KpSC assemblies for ground truth generation.

Strategy:
  1. Download AllTheBacteria v0.2 metadata (maps run_accession → assembly ftp URL).
  2. For each run accession in kpsc_sra_accessions.txt, look up the assembly URL.
  3. Download and decompress the assembly FASTA → data/assemblies/<sample_id>.fasta

AllTheBacteria metadata is a ~1 GB TSV (gzipped).  It is cached locally after
first download.

Usage:
    python data/setup/download_assemblies.py \\
        --accessions   assets/kpsc_sra_accessions.txt \\
        --out-dir      data/assemblies/ \\
        --workers      8

    # Dry-run — just print URLs without downloading:
    python data/setup/download_assemblies.py \\
        --accessions   assets/kpsc_sra_accessions.txt \\
        --out-dir      data/assemblies/ \\
        --dry-run
"""

import argparse
import gzip
import io
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed


ATB_METADATA_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/"
    "Releases/0.2/metadata/allthebacteria_v02_metadata.tsv.gz"
)
ATB_METADATA_CACHE = "assets/atb_metadata_v02.tsv.gz"

# Fallback: ENA portal API for assemblies not in ATB metadata
ENA_FILEREPORT_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run&fields=run_accession,"
    "sample_accession&format=tsv"
)
ENA_ASSEMBLY_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={sample}&result=assembly&fields=accession,"
    "fasta_ftp&format=tsv"
)


# ---------------------------------------------------------------------------
# Metadata loading
# ---------------------------------------------------------------------------

def _download_file(url: str, dest: str, label: str = "") -> None:
    print(f"Downloading {label or url} …", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def load_atb_metadata(cache_path: str) -> dict[str, str]:
    """Return {run_accession: fasta_ftp_url} from AllTheBacteria metadata."""
    if not os.path.exists(cache_path):
        _download_file(ATB_METADATA_URL, cache_path, "AllTheBacteria v0.2 metadata (~1 GB)")

    print(f"Parsing ATB metadata ({cache_path}) …", flush=True)
    run_to_ftp: dict[str, str] = {}

    with gzip.open(cache_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        # Look for run_accession and fasta_ftp (or assembly_ftp) columns
        run_col  = next((i for i, h in enumerate(header) if "run" in h.lower()), None)
        ftp_col  = next(
            (i for i, h in enumerate(header)
             if "fasta" in h.lower() or "assembly" in h.lower()),
            None,
        )
        if run_col is None or ftp_col is None:
            raise RuntimeError(
                f"Cannot find run/fasta columns in ATB metadata. Header: {header}"
            )
        print(f"  run column: '{header[run_col]}'  fasta column: '{header[ftp_col]}'", flush=True)

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(run_col, ftp_col):
                continue
            run = parts[run_col].strip()
            ftp = parts[ftp_col].strip()
            if run and ftp and ftp != "-":
                run_to_ftp[run] = ftp

    print(f"Loaded {len(run_to_ftp):,} run→FASTA mappings from ATB metadata.", flush=True)
    return run_to_ftp


# ---------------------------------------------------------------------------
# ENA API fallback
# ---------------------------------------------------------------------------

def _ena_fasta_url(run_acc: str) -> str | None:
    """Look up assembly FASTA URL via ENA portal API.  Returns None if not found."""
    try:
        with urllib.request.urlopen(ENA_FILEREPORT_URL.format(acc=run_acc), timeout=30) as r:
            lines = r.read().decode().strip().splitlines()
        if len(lines) < 2:
            return None
        parts = lines[1].split("\t")
        if len(parts) < 2:
            return None
        sample_acc = parts[1].strip()
        if not sample_acc:
            return None

        with urllib.request.urlopen(ENA_ASSEMBLY_URL.format(sample=sample_acc), timeout=30) as r:
            lines = r.read().decode().strip().splitlines()
        if len(lines) < 2:
            return None
        ftp = lines[1].split("\t")[-1].strip()
        return ftp if ftp and ftp != "-" else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Download worker
# ---------------------------------------------------------------------------

def _download_one(run_acc: str, ftp_url: str, out_dir: str) -> tuple[str, str]:
    """Download one assembly FASTA.  Returns (run_acc, status)."""
    dest = os.path.join(out_dir, f"{run_acc}.fasta")
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        return run_acc, "skipped"

    # ENA FTP URLs sometimes start with "ftp://" or are bare paths
    url = ftp_url if ftp_url.startswith(("http", "ftp")) else f"ftp://{ftp_url}"

    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            raw = r.read()

        # Decompress if gzipped
        if url.endswith(".gz") or raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)

        with open(dest, "wb") as f:
            f.write(raw)
        return run_acc, "ok"
    except Exception as e:
        return run_acc, f"FAILED: {e}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--accessions", default="assets/kpsc_sra_accessions.txt",
                        help="One SRA run accession per line.")
    parser.add_argument("--out-dir",    default="data/assemblies/",
                        help="Directory for downloaded FASTA files.")
    parser.add_argument("--metadata-cache", default=ATB_METADATA_CACHE,
                        help="Local cache path for ATB metadata TSV.gz.")
    parser.add_argument("--workers",    type=int, default=8,
                        help="Parallel download workers.")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Print URLs without downloading.")
    args = parser.parse_args()

    with open(args.accessions) as f:
        accessions = [l.strip() for l in f if l.strip()]
    print(f"Accessions to process: {len(accessions)}", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load ATB metadata ──────────────────────────────────────────────────
    run_to_ftp = load_atb_metadata(args.metadata_cache)

    # ── Resolve URLs ───────────────────────────────────────────────────────
    tasks: list[tuple[str, str]] = []
    missing: list[str] = []
    for acc in accessions:
        if acc in run_to_ftp:
            tasks.append((acc, run_to_ftp[acc]))
        else:
            missing.append(acc)

    print(f"Found in ATB metadata: {len(tasks)}  |  not found: {len(missing)}", flush=True)

    # ENA API fallback for missing accessions
    if missing:
        print(f"Trying ENA API for {len(missing)} missing accessions …", flush=True)
        for acc in missing:
            url = _ena_fasta_url(acc)
            if url:
                tasks.append((acc, url))
                print(f"  ENA found: {acc} → {url[:60]}…", flush=True)
            else:
                print(f"  WARNING: no assembly found for {acc}", file=sys.stderr)

    print(f"\nTotal assemblies to download: {len(tasks)}", flush=True)

    if args.dry_run:
        for acc, url in tasks:
            print(f"{acc}\t{url}")
        return

    # ── Download ───────────────────────────────────────────────────────────
    n_ok, n_skip, n_fail = 0, 0, 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_download_one, acc, url, args.out_dir): acc
                   for acc, url in tasks}
        for i, fut in enumerate(as_completed(futures), 1):
            acc, status = fut.result()
            if status == "ok":
                n_ok += 1
            elif status == "skipped":
                n_skip += 1
            else:
                n_fail += 1
                print(f"  {status} [{acc}]", flush=True)
            if i % 50 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)} | ok={n_ok} skip={n_skip} fail={n_fail}", flush=True)

    print(f"\nDone. ok={n_ok}  skipped={n_skip}  failed={n_fail}", flush=True)
    print(f"Assemblies in {args.out_dir}: {len(os.listdir(args.out_dir))}", flush=True)


if __name__ == "__main__":
    main()
