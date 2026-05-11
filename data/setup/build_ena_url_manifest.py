"""
build_ena_url_manifest.py
─────────────────────────
For each KpSC expansion accession that does not yet have a BAM, query the
ENA Portal API to obtain the canonical HTTPS FASTQ URLs.

Produces:
    assets/ena_url_manifest.tsv
      accession | layout | r1_url | r2_url
      - layout  : PAIRED or SINGLE
      - r1_url  : HTTPS URL for read 1 (or the single FASTQ for SINGLE)
      - r2_url  : HTTPS URL for read 2 (empty for SINGLE layout)

Usage (run on the HPC login node or locally):
    python3 data/setup/build_ena_url_manifest.py

Requirements: requests (pip install requests)

The script is idempotent: accessions already in the manifest are skipped
on re-runs.
"""

import os, sys, time, csv
import requests

REPO_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ACCESSIONS_FILE = os.path.join(REPO_DIR, "assets", "kpsc_expansion_sra_accessions.txt")
BAM_DIR         = os.path.join(REPO_DIR, "data", "raw", "bam")
OUT_MANIFEST    = os.path.join(REPO_DIR, "assets", "ena_url_manifest.tsv")

BATCH_SIZE  = 200   # ENA API handles up to ~200 per request comfortably
SLEEP_SEC   = 1.0   # polite delay between batches
MAX_RETRIES = 5

ENA_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?result=read_run"
    "&fields=fastq_ftp,library_layout"
    "&format=tsv"
    "&accession={accs}"
)


# ── helpers ───────────────────────────────────────────────────────────────────

def ftp_to_https(ftp_path: str) -> str:
    """Convert bare FTP path → HTTPS URL (add https:// prefix)."""
    if not ftp_path:
        return ""
    path = ftp_path.strip()
    if path.startswith("ftp://") or path.startswith("https://"):
        return path.replace("ftp://", "https://")
    return "https://" + path


def query_ena(accessions: list[str]) -> dict[str, dict]:
    """
    Query ENA Portal API for a batch of accessions.
    Returns dict: accession → {layout, r1_url, r2_url}
    """
    url = ENA_API.format(accs=",".join(accessions))
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, timeout=60)
            if resp.status_code == 200:
                break
            print(f"  HTTP {resp.status_code} on attempt {attempt}; retrying…")
        except requests.RequestException as e:
            print(f"  Request error on attempt {attempt}: {e}; retrying…")
        time.sleep(attempt * 5)
    else:
        print(f"  ERROR: all {MAX_RETRIES} attempts failed for batch — skipping.")
        return {}

    results = {}
    reader = csv.DictReader(resp.text.splitlines(), delimiter="\t")
    for row in reader:
        acc = row.get("run_accession", "").strip()
        if not acc:
            continue
        ftp_field = row.get("fastq_ftp", "").strip()
        layout    = (row.get("library_layout") or "PAIRED").strip().upper()

        urls = [ftp_to_https(u) for u in ftp_field.split(";") if u.strip()]

        # Separate R1 / R2 / single
        r1 = r2 = ""
        single = ""
        for u in urls:
            fname = u.rsplit("/", 1)[-1]
            if fname.endswith("_1.fastq.gz"):
                r1 = u
            elif fname.endswith("_2.fastq.gz"):
                r2 = u
            elif fname.endswith(".fastq.gz") and "_" not in fname.rsplit("/", 1)[-1].replace(".fastq.gz", ""):
                single = u

        if layout == "PAIRED" and r1 and r2:
            results[acc] = {"layout": "PAIRED", "r1_url": r1, "r2_url": r2}
        elif r1 and r2:
            # API says SINGLE but we have paired files — trust the files
            results[acc] = {"layout": "PAIRED", "r1_url": r1, "r2_url": r2}
        elif single:
            results[acc] = {"layout": "SINGLE", "r1_url": single, "r2_url": ""}
        elif r1:
            results[acc] = {"layout": "SINGLE", "r1_url": r1, "r2_url": ""}
        else:
            print(f"  WARNING: no FASTQ URLs found for {acc} (ftp field: {ftp_field!r})")

    return results


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    # Load all accessions
    with open(ACCESSIONS_FILE) as fh:
        all_accs = [line.strip() for line in fh if line.strip()]
    print(f"Total accessions: {len(all_accs):,}")

    # Which ones already have BAMs?
    done_bams = {
        os.path.basename(f).replace(".bam", "")
        for f in os.listdir(BAM_DIR)
        if f.endswith(".bam")
    } if os.path.isdir(BAM_DIR) else set()
    print(f"Already have BAMs: {len(done_bams):,}")

    # Which ones are already in the manifest?
    existing_manifest: dict[str, dict] = {}
    if os.path.exists(OUT_MANIFEST):
        with open(OUT_MANIFEST) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                existing_manifest[row["accession"]] = row
        print(f"Already in manifest: {len(existing_manifest):,}")

    # Need URL for: missing BAM + not yet in manifest
    need_url = [
        a for a in all_accs
        if a not in done_bams and a not in existing_manifest
    ]
    print(f"Need ENA API query: {len(need_url):,}")

    if not need_url:
        print("Nothing to do — manifest is complete.")
        return

    # Query ENA in batches
    new_results: dict[str, dict] = {}
    total_batches = (len(need_url) + BATCH_SIZE - 1) // BATCH_SIZE
    for i in range(0, len(need_url), BATCH_SIZE):
        batch = need_url[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        print(f"  Batch {batch_num}/{total_batches}: {batch[0]} … {batch[-1]}", end="", flush=True)
        result = query_ena(batch)
        new_results.update(result)
        print(f"  → {len(result)}/{len(batch)} resolved")
        time.sleep(SLEEP_SEC)

    # Write manifest (append new rows; rewrite if file exists for header consistency)
    all_results = {**existing_manifest, **new_results}
    with open(OUT_MANIFEST, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["accession", "layout", "r1_url", "r2_url"],
            delimiter="\t",
        )
        writer.writeheader()
        for acc in all_accs:
            if acc in all_results:
                row = all_results[acc]
                writer.writerow({
                    "accession": acc,
                    "layout":    row.get("layout", ""),
                    "r1_url":    row.get("r1_url", ""),
                    "r2_url":    row.get("r2_url", ""),
                })

    # Summary
    paired  = sum(1 for r in all_results.values() if r.get("layout") == "PAIRED")
    single  = sum(1 for r in all_results.values() if r.get("layout") == "SINGLE")
    missing = len(need_url) - len(new_results)

    print(f"\nManifest written: {OUT_MANIFEST}")
    print(f"  PAIRED  : {paired:,}")
    print(f"  SINGLE  : {single:,}")
    if missing:
        print(f"  NOT FOUND in ENA: {missing:,}  (accessions with no FASTQ listed)")


if __name__ == "__main__":
    main()
