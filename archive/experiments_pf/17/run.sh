#!/bin/bash
# 04_gene_cnv_caller adds a CRR gate for long genes (CRT, MDR1, PM2_PM3); HMM unchanged.
# Copies all outputs from exp 16, then re-runs CNV calling and evaluation only.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

EXP16="$REPO_ROOT/data/results/16_skip_st0.80_crrthresh1.40"
OUT_DIR="$REPO_ROOT/data/results/17_skip_st0.80_crr_gate"

mkdir -p "$OUT_DIR"
cp "$EXP16/reconstructions.npy" "$OUT_DIR/"
cp "$EXP16/sample_ids.npy"      "$OUT_DIR/"
cp "$EXP16/latents.npy"         "$OUT_DIR/"
cp "$EXP16/segments.parquet"    "$OUT_DIR/"

"$REPO_ROOT/.venv/bin/python" - <<'PYEOF'
import sys, os, yaml, importlib
sys.path.insert(0, os.getcwd())

config_path = "experiments/17/config.yaml"
config_dir  = os.path.dirname(os.path.abspath(config_path))

with open(config_path) as f:
    cfg = yaml.safe_load(f)

def resolve(path):
    return path if os.path.isabs(path) else os.path.join(config_dir, path)

store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_cnv_calls  = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

print("Calling gene CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg)

cfg_resolved = dict(cfg)
cfg_resolved["pf9_gt_path"]   = resolve(cfg["pf9_gt_path"])
cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
print("Running evaluation...", flush=True)
run_evaluation(out_dir, cfg_resolved)
PYEOF
