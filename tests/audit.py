#!/usr/bin/env python3
"""
CNVRock pipeline audit — quick sanity checks across environments,
reference files, data integrity, and results consistency.

Run locally for local checks:
    python3 tests/audit.py

Run on HPC for full check:
    /home/lshlt19/miniconda3/envs/cnvrock/bin/python3 tests/audit.py --hpc
"""
from __future__ import annotations
import argparse, ast, shutil, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
HPC_REPO = Path("/home/lshlt19/CNVRock")

PASS = "\033[32m✓\033[0m"
FAIL = "\033[31m✗\033[0m"
WARN = "\033[33m~\033[0m"
SKIP = "\033[90m–\033[0m"

results: list[tuple[str, str, str]] = []   # (status, category, message)


def ok(cat, msg):
    results.append(("PASS", cat, msg))
    print(f"  {PASS} [{cat}] {msg}")


def fail(cat, msg):
    results.append(("FAIL", cat, msg))
    print(f"  {FAIL} [{cat}] {msg}")


def warn(cat, msg):
    results.append(("WARN", cat, msg))
    print(f"  {WARN} [{cat}] {msg}")


def skip(cat, msg):
    results.append(("SKIP", cat, msg))
    print(f"  {SKIP} [{cat}] {msg}")


# ── Helpers ──────────────────────────────────────────────────────────────────

def check_file(path: Path, cat: str, min_bytes: int = 1,
               min_lines: int = 0) -> bool:
    if not path.exists():
        fail(cat, f"MISSING: {path.name}")
        return False
    size = path.stat().st_size
    if size < min_bytes:
        fail(cat, f"EMPTY/TOO_SMALL: {path.name} ({size} bytes)")
        return False
    if min_lines:
        try:
            n = sum(1 for _ in open(path))
            if n < min_lines:
                fail(cat, f"TOO_FEW_LINES: {path.name} ({n} < {min_lines})")
                return False
        except Exception:
            pass
    ok(cat, f"{path.name} exists ({size/1024:.0f} KB)")
    return True


def check_tsv_columns(path: Path, cat: str,
                      required_cols: list[str]) -> bool:
    if not path.exists():
        fail(cat, f"MISSING: {path.name}")
        return False
    try:
        header = path.open().readline().rstrip().split("\t")
        missing = [c for c in required_cols if c not in header]
        if missing:
            fail(cat, f"{path.name}: missing columns {missing}")
            return False
        ok(cat, f"{path.name}: columns OK {required_cols}")
        return True
    except Exception as e:
        fail(cat, f"{path.name}: {e}")
        return False


def check_fasta_contigs(path: Path, cat: str,
                        expected: list[str]) -> None:
    if not path.exists():
        fail(cat, f"MISSING: {path.name}"); return
    found = set()
    for line in open(path):
        if line.startswith(">"):
            found.add(line[1:].split()[0].strip())
    missing = [c for c in expected if c not in found]
    if missing:
        fail(cat, f"{path.name}: missing contigs {missing}")
    else:
        ok(cat, f"{path.name}: all {len(expected)} expected contigs present")


def check_metric(eval_path: Path, cat: str,
                 gene: str, metric: str, expected: float,
                 tol: float = 0.02) -> None:
    """
    Parse metric from evaluation.txt table rows like:
        blaKPC       presence   0.98  0.04  1.00       1.00     6261
    Columns: gene  mode  MCC  FNR  PPV  call_rate  n_eval
    """
    if not eval_path.exists():
        skip(cat, f"evaluation.txt not found — {gene} {metric} not checked")
        return
    COL = {"MCC": 2, "FNR": 3, "PPV": 4, "call_rate": 5}
    col_idx = COL.get(metric.upper(), 2)
    for line in open(eval_path):
        parts = line.split()
        if parts and parts[0] == gene and len(parts) >= col_idx + 1:
            try:
                val = float(parts[col_idx])
            except ValueError:
                continue
            if abs(val - expected) <= tol:
                ok(cat, f"{gene} {metric}={val:.3f} (expected ~{expected})")
            else:
                warn(cat, f"{gene} {metric}={val:.3f} "
                          f"(expected ~{expected}, diff={abs(val-expected):.3f})")
            return
    skip(cat, f"{gene} {metric} not found in {eval_path.name}")


# ── Test suites ───────────────────────────────────────────────────────────────

def test_syntax(repo: Path) -> None:
    """All Python scripts parse without syntax errors."""
    print("\n── Syntax (Python AST parse) ────────────────────────────────")
    py_dirs = ["analysis", "models", "data/setup", "tests"]
    errors = []
    n_ok = 0
    for d in py_dirs:
        for f in sorted((repo / d).rglob("*.py")):
            try:
                ast.parse(f.read_text())
                n_ok += 1
            except SyntaxError as e:
                errors.append((f.relative_to(repo), e))
    if errors:
        for f, e in errors:
            fail("syntax", f"{f}: {e}")
    else:
        ok("syntax", f"All {n_ok} Python files parse cleanly")


def test_tools_local(repo: Path) -> None:
    """Tools available locally for building manuscript."""
    print("\n── Local tools ──────────────────────────────────────────────")
    is_hpc = str(repo) == str(HPC_REPO)
    # python-docx — only needed locally (manuscript.docx built on laptop)
    try:
        import docx  # noqa
        ok("tools-local", "python-docx importable")
    except ImportError:
        if is_hpc:
            skip("tools-local", "python-docx not on HPC (only needed locally for docx build)")
        else:
            fail("tools-local", "python-docx NOT importable — manuscript.docx cannot be built")
    # pyarrow (for ATB parquet)
    try:
        import pyarrow  # noqa
        ok("tools-local", "pyarrow importable (ATB parquet access)")
    except ImportError:
        warn("tools-local", "pyarrow not available — ATB queries will fail locally")


def test_tools_hpc() -> None:
    """Tools available in expected HPC environments."""
    print("\n── HPC tool paths ───────────────────────────────────────────")
    checks = [
        ("/home/lshlt19/miniconda3/envs/cnvrock/bin/python3",  "cnvrock Python"),
        ("/home/lshlt19/miniconda3/envs/aligners/bin/samtools","samtools (aligners)"),
        ("/home/lshlt19/miniconda3/envs/blast_env/bin/blastn", "blastn (blast_env)"),
        ("/home/lshlt19/miniconda3/envs/blast_env/bin/makeblastdb","makeblastdb (blast_env)"),
    ]
    for path_str, label in checks:
        p = Path(path_str)
        if p.exists():
            ok("tools-hpc", f"{label}: {path_str}")
        else:
            # try shutil.which as fallback
            found = shutil.which(p.name)
            if found:
                warn("tools-hpc", f"{label}: not at expected path, but found at {found}")
            else:
                fail("tools-hpc", f"{label}: NOT FOUND at {path_str}")


def test_reference_files(repo: Path) -> None:
    """Key reference FASTA and coordinate files.
    Large FASTAs live only on HPC — skip locally."""
    print("\n── Reference files ──────────────────────────────────────────")
    is_hpc = str(repo) == str(HPC_REPO)
    # KpSC
    if is_hpc:
        check_fasta_contigs(
            repo / "assets/HS11286_extended.fasta", "ref-kpsc",
            ["NC_016845.1", "NC_016846.1", "MK552109.1", "MZ606384.2",
             "JN626286.1", "MH287085.1", "MN540571.1", "MZ382871.1"])
    else:
        skip("ref-kpsc", "HS11286_extended.fasta (HPC-only, gitignored)")

    coords = repo / "assets/plasmid_refs/plasmid_gene_coords.tsv"
    if is_hpc:
        check_tsv_columns(coords, "ref-kpsc", ["gene", "contig", "start", "end"])
        if coords.exists():
            txt = coords.read_text()
            if "aac3-II" in txt:
                fail("ref-kpsc", "plasmid_gene_coords.tsv still has aac3-II (should be removed)")
            else:
                ok("ref-kpsc", "aac3-II absent from plasmid_gene_coords.tsv ✓")
    else:
        skip("ref-kpsc", "plasmid_gene_coords.tsv (HPC-only)")

    # A. baumannii
    if is_hpc:
        check_fasta_contigs(
            repo / "assets/abaumannii_ref/AB5075.fasta", "ref-abaum",
            ["NZ_CP008706.1"])
    else:
        skip("ref-abaum", "AB5075.fasta (HPC-only)")
    check_tsv_columns(
        repo / "assets/abaumannii_ref/gene_coords.tsv", "ref-abaum",
        ["gene", "contig", "start", "end"])

    # C. auris
    if is_hpc:
        check_fasta_contigs(
            repo / "assets/cauris_ref/B8441v3.fasta", "ref-cauris",
            ["CM076438.1", "CM076439.1", "CM076440.1", "CM076442.1"])
    else:
        skip("ref-cauris", "B8441v3.fasta (HPC-only)")
    check_tsv_columns(
        repo / "assets/cauris_ref/gene_coords.tsv", "ref-cauris",
        ["gene", "contig", "start", "end"])


def test_ground_truth_files(repo: Path) -> None:
    """AMRFinder GT files exist and have expected structure."""
    print("\n── Ground-truth files ───────────────────────────────────────")
    is_hpc = str(repo) == str(HPC_REPO)

    check_tsv_columns(
        repo / "assets/amrfinder_gt_expansion.tsv", "gt",
        ["sample_id", "biosample", "blaKPC", "blaNDM", "blaCTX-M"])

    pe = repo / "assets/amrfinder_gt_expansion_phaseE.tsv"
    if is_hpc:
        check_tsv_columns(pe, "gt", ["sample_id", "sul1", "sul2", "dfrA12", "dfrA14"])
        if pe.exists():
            header = pe.open().readline()
            if "aac3-II" in header:
                fail("gt", "amrfinder_gt_expansion_phaseE.tsv has aac3-II column (should be removed)")
            else:
                ok("gt", "aac3-II correctly absent from phaseE GT ✓")
    else:
        skip("gt", "amrfinder_gt_expansion_phaseE.tsv (HPC-only)")

    check_file(
        repo / "assets/abaumannii_amrfinder_gt.tsv", "gt",
        min_lines=1450)
    check_tsv_columns(
        repo / "assets/abaumannii_amrfinder_gt.tsv", "gt",
        ["sample_id", "blaOXA-23", "blaOXA-24-like", "blaOXA-58-like"])


def test_experiment_configs(repo: Path) -> None:
    """Experiment configs point to existing paths."""
    print("\n── Experiment configs ───────────────────────────────────────")
    import yaml  # may not be available
    exp_dirs = sorted((repo / "models/experiments").glob("*/"))
    for exp_dir in exp_dirs:
        cfg_path = exp_dir / "config.yaml"
        if not cfg_path.exists():
            continue
        try:
            cfg = yaml.safe_load(cfg_path.read_text())
        except Exception as e:
            warn("exp-cfg", f"{exp_dir.name}: YAML parse error: {e}")
            continue
        exp_name = exp_dir.name
        # Check chrom_gene_coords_path if set
        gcpath = cfg.get("chrom_gene_coords_path")
        if gcpath:
            resolved = (exp_dir / gcpath).resolve()
            if not resolved.exists():
                fail("exp-cfg", f"{exp_name}: chrom_gene_coords_path not found: {resolved}")
            else:
                ok("exp-cfg", f"{exp_name}: chrom_gene_coords_path OK")
        # Check evaluation module exists
        eval_mod = cfg.get("evaluation")
        if eval_mod:
            eval_file = repo / "models/evaluation" / f"{eval_mod}.py"
            if not eval_file.exists():
                fail("exp-cfg", f"{exp_name}: evaluation module missing: {eval_mod}.py")
            else:
                ok("exp-cfg", f"{exp_name}: evaluation={eval_mod} OK")


def test_results_consistency(repo: Path) -> None:
    """Check headline metrics in evaluation files match manuscript."""
    print("\n── Results consistency (manuscript numbers) ─────────────────")
    is_hpc = str(repo) == str(HPC_REPO)
    eval_33 = repo / "data/results/33_kpsc_expansion_10k/evaluation.txt"
    check_metric(eval_33, "results", "blaKPC",   "MCC",       0.98)
    check_metric(eval_33, "results", "blaNDM",   "MCC",       0.76)
    check_metric(eval_33, "results", "blaCTX-M", "MCC",       0.67)
    check_metric(eval_33, "results", "aac6-Ib-cr","MCC",      0.84)

    eval_37 = repo / "data/results/37_kpsc_phase_e_10k/evaluation.txt"
    check_metric(eval_37, "results", "dfrA14",   "MCC",       0.88)
    check_metric(eval_37, "results", "dfrA12",   "MCC",       0.81)
    check_metric(eval_37, "results", "sul2",     "MCC",       0.81)

    # exp40 — Phase 1 hold-out OOD evaluation (§3.7)
    eval_40 = repo / "data/results/40_phase1_holdout_10k/evaluation.txt"
    check_metric(eval_40, "results", "blaKPC",   "MCC",       1.00, tol=0.01)
    check_metric(eval_40, "results", "blaNDM",   "MCC",       0.99, tol=0.02)
    check_metric(eval_40, "results", "blaCTX-M", "MCC",       0.82, tol=0.03)
    check_metric(eval_40, "results", "aac6-Ib-cr","MCC",      0.86, tol=0.03)

    # C. auris mutation GT (HPC path)
    mut_gt = repo / "data/results/cauris_mutation_gt/cauris_erg11_mutation_gt.tsv"
    if not is_hpc:
        skip("results", "cauris ERG11 GT (HPC-only)")
    elif mut_gt.exists():
        try:
            import pandas as pd
        except ImportError as e:
            warn("results", f"pandas unavailable (GLIBCXX mismatch?): {e}")
            return
        df = pd.read_csv(mut_gt, sep="\t")
        n_mut = int(df["erg11_R_mutation"].sum()) if "erg11_R_mutation" in df.columns else 0
        if n_mut > 100:
            ok("results", f"C. auris ERG11 mutations: {n_mut}/503 (>100, plausible)")
        elif n_mut == 0:
            fail("results", "C. auris ERG11 mutations: 0 — samtools path bug? "
                            "Run with /home/lshlt19/miniconda3/envs/aligners/bin/samtools")
        else:
            warn("results", f"C. auris ERG11 mutations: {n_mut} (low — check)")
    else:
        fail("results", "cauris_erg11_mutation_gt.tsv not found on HPC")

    # VAE ablation
    abl = repo / "data/results/vae_ablation/ablation_summary.tsv"
    if abl.exists():
        try:
            import pandas as pd
        except ImportError:
            skip("results", "VAE ablation check skipped (pandas unavailable)")
            return
        df = pd.read_csv(abl, sep="\t")
        vae = df[df["baseline"] == "C_vae"]
        if not vae.empty:
            rmse = float(vae["spikein_recovery_RMSE"].iloc[0])
            if rmse < 0.75:
                ok("results", f"VAE spike-in RMSE={rmse:.4f} (best baseline)")
            else:
                warn("results", f"VAE spike-in RMSE={rmse:.4f} (unexpectedly high)")
    else:
        skip("results", "VAE ablation summary not found")


def test_manuscript(repo: Path) -> None:
    """Manuscript sanity checks."""
    print("\n── Manuscript ───────────────────────────────────────────────")
    md = repo / "paper/manuscript.md"
    if not md.exists():
        skip("manuscript", "manuscript.md not found (gitignored — check locally)")
        return
    text = md.read_text()
    # No VERIFY tags
    n_verify = text.count("[VERIFY")
    if n_verify == 0:
        ok("manuscript", "No [VERIFY] tags remaining")
    else:
        fail("manuscript", f"{n_verify} [VERIFY] tags still present")
    # No placeholders
    n_ph = text.count("⟨")
    if n_ph <= 1:   # §3.7 hold-out is the only known remaining placeholder
        ok("manuscript", f"Only {n_ph} placeholder(s) remaining (§3.7 hold-out expected)")
    else:
        warn("manuscript", f"{n_ph} placeholders remaining")
    # References exist
    import re
    refs = set(re.findall(r"\[(\d+)\]", text))
    if refs:
        max_ref = max(int(r) for r in refs)
        ok("manuscript", f"References [1]–[{max_ref}] cited in text")
    # Title check
    if "cross-kingdom" in text.lower():
        ok("manuscript", "Title contains 'cross-kingdom'")
    else:
        warn("manuscript", "Title may not be updated")
    # .docx exists
    check_file(repo / "paper/manuscript.docx", "manuscript", min_bytes=50000)
    check_file(repo / "paper/psi_nature_microbiology.docx", "manuscript", min_bytes=20000)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="CNVRock audit")
    parser.add_argument("--hpc", action="store_true",
                        help="Run on HPC (use /home/lshlt19/CNVRock paths)")
    args = parser.parse_args()

    repo = HPC_REPO if args.hpc else REPO
    print(f"\n{'='*60}")
    print(f"CNVRock Audit  |  repo={repo}  |  {'HPC' if args.hpc else 'local'}")
    print(f"{'='*60}")

    test_syntax(REPO)          # always uses local repo for syntax check
    test_tools_local(repo)

    if args.hpc:
        test_tools_hpc()
    else:
        skip("tools-hpc", "Skipped (run with --hpc on HPC)")

    test_reference_files(repo)
    test_ground_truth_files(repo)

    try:
        test_experiment_configs(repo)
    except ImportError:
        skip("exp-cfg", "PyYAML not available — skipping config checks")

    test_results_consistency(repo)
    test_manuscript(repo)

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    n_pass = sum(1 for s, *_ in results if s == "PASS")
    n_fail = sum(1 for s, *_ in results if s == "FAIL")
    n_warn = sum(1 for s, *_ in results if s == "WARN")
    n_skip = sum(1 for s, *_ in results if s == "SKIP")
    print(f"PASS={n_pass}  FAIL={n_fail}  WARN={n_warn}  SKIP={n_skip}")
    if n_fail > 0:
        print(f"\n{FAIL} FAILURES:")
        for s, cat, msg in results:
            if s == "FAIL":
                print(f"   [{cat}] {msg}")
        sys.exit(1)
    else:
        print(f"\n{PASS} All critical checks passed.")


if __name__ == "__main__":
    main()
