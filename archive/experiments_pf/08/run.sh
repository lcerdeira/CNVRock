#!/bin/bash
# Evaluation-only: copy exp 07 gene_calls.tsv and segments.parquet into a new
# output directory, then run 03_pf9_evaluation to add segment diagnostics.
# No training, HMM re-segmentation, or CNV calling.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
PYTHON="$REPO_ROOT/.venv/bin/python"
SRC="$REPO_ROOT/data/results/07_higher_self_transition"
OUT="$REPO_ROOT/data/results/08_new_evaluation"

mkdir -p "$OUT"
cp "$SRC/gene_calls.tsv"   "$OUT/"
cp "$SRC/segments.parquet" "$OUT/"

REPO_ROOT="$REPO_ROOT" "$PYTHON" - <<'PYEOF'
import os, sys, yaml, importlib

repo = os.environ['REPO_ROOT']
sys.path.insert(0, os.path.join(repo, 'models'))

cfg_path = os.path.join(repo, 'models/experiments/08/config.yaml')
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)

config_dir = os.path.dirname(cfg_path)
def resolve(p):
    return p if os.path.isabs(p) else os.path.join(config_dir, p)

out_dir = resolve(cfg['out_dir'])
cfg['pf9_gt_path']   = resolve(cfg['pf9_gt_path'])
if cfg.get('pf9_meta_path'):
    cfg['pf9_meta_path'] = resolve(cfg['pf9_meta_path'])

run_evaluation = importlib.import_module('evaluation.' + cfg['evaluation']).run_evaluation
run_evaluation(out_dir, cfg)
PYEOF
