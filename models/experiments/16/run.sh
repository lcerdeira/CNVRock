#!/bin/bash
# Only cnv_crr_amp_threshold changed (1.30 → 1.40); HMM parameters and architecture unchanged.
# Copies all outputs from exp 15 (including segments.parquet), then re-runs CNV calling
# and evaluation only. HMM step skipped — parameters are identical to exp 15.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

EXP15="$REPO_ROOT/data/results/15_skip_connections_st0.80"
OUT_DIR="$REPO_ROOT/data/results/16_skip_st0.80_crrthresh1.40"

mkdir -p "$OUT_DIR"
cp "$EXP15/reconstructions.npy" "$OUT_DIR/"
cp "$EXP15/sample_ids.npy"      "$OUT_DIR/"
cp "$EXP15/latents.npy"         "$OUT_DIR/"
cp "$EXP15/segments.parquet"    "$OUT_DIR/"

"$REPO_ROOT/.venv/bin/python" - <<'PYEOF'
import sys, os, yaml, importlib
sys.path.insert(0, os.getcwd())

config_path = "experiments/16/config.yaml"
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
