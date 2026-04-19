#!/bin/bash
# Reuse exp 14 (05_conv_vae) inference outputs; only HMM self_transition changed (0.65 → 0.80).
# Copies reconstructions/sample_ids/latents from exp 14, then re-runs HMM, CNV calling,
# and evaluation with the new config. No retraining or re-inference required.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$REPO_ROOT/models"

EXP14="$REPO_ROOT/data/results/14_skip_connections"
OUT_DIR="$REPO_ROOT/data/results/15_skip_connections_st0.80"

mkdir -p "$OUT_DIR"
cp "$EXP14/reconstructions.npy" "$OUT_DIR/"
cp "$EXP14/sample_ids.npy"      "$OUT_DIR/"
cp "$EXP14/latents.npy"         "$OUT_DIR/"

"$REPO_ROOT/.venv/bin/python" - <<'PYEOF'
import sys, os, yaml, importlib
sys.path.insert(0, os.getcwd())

config_path = "experiments/15/config.yaml"
config_dir  = os.path.dirname(os.path.abspath(config_path))

with open(config_path) as f:
    cfg = yaml.safe_load(f)

def resolve(path):
    return path if os.path.isabs(path) else os.path.join(config_dir, path)

store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_hmm_all_samples = importlib.import_module(f"hmm.{cfg['hmm']}").run_hmm_all_samples
run_cnv_calls       = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation      = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

print("Fitting HMM segments...", flush=True)
run_hmm_all_samples(store_path, out_dir, cfg)

print("Calling gene CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg)

cfg_resolved = dict(cfg)
cfg_resolved["pf9_gt_path"]  = resolve(cfg["pf9_gt_path"])
cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
print("Running evaluation...", flush=True)
run_evaluation(out_dir, cfg_resolved)
PYEOF
