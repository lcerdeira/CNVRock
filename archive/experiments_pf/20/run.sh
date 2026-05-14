#!/bin/bash
# Exp 20: lower GCH1 fallback threshold on exp 19 CNV caller.
# Symlinks exp 19 inference, HMM, and VAE outputs into the new out_dir, then
# runs only CNV calling + evaluation (no retraining, no re-segmentation).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP19_OUT="$REPO_ROOT/data/results/19_skip_st0.80_sin_partialspan"
EXP20_OUT="$REPO_ROOT/data/results/20_skip_st0.80_sin_partialspan_gate1.25"

mkdir -p "$EXP20_OUT"
ln -sf "$EXP19_OUT/checkpoint.pth"       "$EXP20_OUT/checkpoint.pth"
ln -sf "$EXP19_OUT/latents.npy"          "$EXP20_OUT/latents.npy"
ln -sf "$EXP19_OUT/reconstructions.npy"  "$EXP20_OUT/reconstructions.npy"
ln -sf "$EXP19_OUT/sample_ids.npy"       "$EXP20_OUT/sample_ids.npy"
ln -sf "$EXP19_OUT/segments.parquet"     "$EXP20_OUT/segments.parquet"

cd "$REPO_ROOT/models"
"$REPO_ROOT/.venv/bin/python" - <<'EOF'
import sys, os, yaml, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

config_path = "experiments/20/config.yaml"
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
