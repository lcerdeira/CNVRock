#!/bin/bash
# Exp 19: fractional gene coverage CNV caller on exp 18 VAE + HMM segments.
# Symlinks exp 18 inference and HMM outputs into the new out_dir, then
# runs only CNV calling + evaluation (no retraining, no re-segmentation).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP18_OUT="$REPO_ROOT/data/results/18_skip_st0.80_sin_finetune"
EXP19_OUT="$REPO_ROOT/data/results/19_skip_st0.80_sin_partialspan"

mkdir -p "$EXP19_OUT"
ln -sf "$EXP18_OUT/checkpoint.pth"       "$EXP19_OUT/checkpoint.pth"
ln -sf "$EXP18_OUT/latents.npy"          "$EXP19_OUT/latents.npy"
ln -sf "$EXP18_OUT/reconstructions.npy"  "$EXP19_OUT/reconstructions.npy"
ln -sf "$EXP18_OUT/sample_ids.npy"       "$EXP19_OUT/sample_ids.npy"
ln -sf "$EXP18_OUT/segments.parquet"     "$EXP19_OUT/segments.parquet"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, os, yaml, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

config_path = "experiments/19/config.yaml"
config_dir  = os.path.dirname(os.path.abspath(config_path))
with open(config_path) as f:
    cfg = yaml.safe_load(f)

def resolve(p):
    return p if os.path.isabs(p) else os.path.join(config_dir, p)

store_path = resolve(cfg["store_path"])
out_dir    = resolve(cfg["out_dir"])

run_cnv_calls  = importlib.import_module(f"cnv.{cfg['cnv']}").run_cnv_calls
run_evaluation = importlib.import_module(f"evaluation.{cfg['evaluation']}").run_evaluation

print("Calling gene CNVs...", flush=True)
run_cnv_calls(store_path, out_dir, cfg)

cfg_resolved = dict(cfg)
cfg_resolved["pf9_gt_path"]  = resolve(cfg["pf9_gt_path"])
cfg_resolved["pf9_meta_path"] = resolve(cfg["pf9_meta_path"])
run_evaluation(out_dir, cfg_resolved)
print("Done.", flush=True)
EOF
