#!/bin/bash
# Exp 29: skip retraining + HMM; re-run plasmid CNV calls (new threshold) + eval.
#
# Changes vs exp 28:
#   - aac6-Ib-cr absent_threshold lowered to 0.10 in plasmid_gene_coords.tsv
#   - GT re-aggregated with expanded qnrB and OXA-48-like patterns
#
# Usage (from repo root):
#   bash models/experiments/29/run_calls_eval.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")/../../.." && pwd)"
EXP27_OUT="$REPO_DIR/data/results/27_kpsc_phase_c_v1"
EXP29_OUT="$REPO_DIR/data/results/29_kpsc_phase_c_v3"
CONFIG="$REPO_DIR/models/experiments/29/config.yaml"

echo "=== Exp 29: plasmid calls + evaluation run ==="
echo "Source (chromosomal outputs): $EXP27_OUT"
echo "Target: $EXP29_OUT"

mkdir -p "$EXP29_OUT"

# Copy chromosomal outputs from exp 27 (no retraining or HMM re-run needed)
for f in checkpoint.pth latents.npy reconstructions.npy sample_ids.npy \
          segments.parquet gene_calls.tsv training_log.json; do
    src="$EXP27_OUT/$f"
    if [[ -f "$src" ]]; then
        cp "$src" "$EXP29_OUT/$f"
        echo "  Copied $f"
    else
        echo "  WARNING: $src not found — skipping"
    fi
done

echo ""
echo "=== Re-running plasmid CNV calls (updated aac6-Ib-cr threshold) ==="
cd "$REPO_DIR"
python - <<'EOF'
import importlib, os, sys, yaml

sys.path.insert(0, os.path.join(os.getcwd(), "models"))

config_path = "models/experiments/29/config.yaml"
config_dir  = os.path.dirname(os.path.abspath(config_path))

with open(config_path) as f:
    cfg = yaml.safe_load(f)

def resolve(p):
    return p if os.path.isabs(p) else os.path.join(config_dir, p)

out_dir      = resolve(cfg["out_dir"])
store_path   = resolve(cfg["store_path"])

run_plasmid_cnv_calls = importlib.import_module(
    f"cnv.{cfg.get('plasmid_cnv', '07_plasmid_cnv_caller')}"
).run_plasmid_cnv_calls

plasmid_cfg = dict(cfg)
plasmid_cfg["store_path"]               = store_path
plasmid_cfg["plasmid_store_path"]       = resolve(cfg["plasmid_store_path"])
plasmid_cfg["plasmid_gene_coords_path"] = resolve(cfg["plasmid_gene_coords_path"])

run_plasmid_cnv_calls(out_dir, plasmid_cfg)
EOF

echo ""
echo "=== Running evaluation ==="
python - <<'EOF'
import importlib, os, sys, yaml

sys.path.insert(0, os.path.join(os.getcwd(), "models"))

config_path = "models/experiments/29/config.yaml"
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
echo "Done. Results in $EXP29_OUT/evaluation.txt"
