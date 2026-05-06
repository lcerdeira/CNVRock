#!/bin/bash
# Exp 28: skip retraining; copy exp 27 model outputs and run evaluation only.
#
# Retraining is unnecessary because the model and plasmid calls are unchanged.
# Only the ground truth TSV is updated (new qnrB1, blaOXA-48, aac6-Ib-cr columns).
#
# Usage (from repo root):
#   bash models/experiments/28/run_eval_only.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP27_OUT="$REPO_DIR/data/results/27_kpsc_phase_c_v1"
EXP28_OUT="$REPO_DIR/data/results/28_kpsc_phase_c_v2"
CONFIG="$REPO_DIR/models/experiments/28/config.yaml"

echo "=== Exp 28: evaluation-only run ==="
echo "Source: $EXP27_OUT"
echo "Target: $EXP28_OUT"

mkdir -p "$EXP28_OUT"

# Copy all inference/call outputs from exp 27 (no training or HMM re-run needed)
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy \
          segments.parquet gene_calls.tsv plasmid_gene_calls.tsv training_log.json; do
    src="$EXP27_OUT/$f"
    if [[ -f "$src" ]]; then
        cp "$src" "$EXP28_OUT/$f"
        echo "  Copied $f"
    else
        echo "  WARNING: $src not found — skipping"
    fi
done

echo ""
echo "=== Running evaluation ==="
cd "$REPO_DIR"
.venv/bin/python - <<'EOF'
import importlib, os, sys, yaml

sys.path.insert(0, os.path.join(os.getcwd(), "models"))

config_path = "models/experiments/28/config.yaml"
config_dir  = os.path.dirname(os.path.abspath(config_path))

with open(config_path) as f:
    cfg = yaml.safe_load(f)

def resolve(p):
    return p if os.path.isabs(p) else os.path.join(config_dir, p)

out_dir = resolve(cfg["out_dir"])

run_evaluation = importlib.import_module(
    f"evaluation.{cfg['evaluation']}"
).run_evaluation

cfg_resolved = dict(cfg)
cfg_resolved["kpsc_gt_path"] = resolve(cfg["kpsc_gt_path"])
for key in ("kpsc_kleborate_gt_path", "kpsc_meta_path", "plasmid_gene_coords_path"):
    if cfg.get(key):
        cfg_resolved[key] = resolve(cfg[key])

run_evaluation(out_dir, cfg_resolved)
EOF

echo ""
echo "Done. Results in $EXP28_OUT/evaluation.txt"
