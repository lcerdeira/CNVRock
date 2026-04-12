# experiments/

Each sub-folder is a self-contained experiment. The `config.yaml` is the single source of truth for everything: data paths, component versions, and all hyperparameters.

```
01/
  config.yaml   data paths, component versions, all hyperparameters
  run.sh        run locally (MPS or CPU)
  experiment.py thin Python wrapper (calls train.py) — prefer run.sh
```

## Running

```bash
cd 01
bash run.sh
```

## Adding a new experiment

```bash
cp -r 01 02
# edit 02/config.yaml — at minimum change out_dir
```

**Reusing existing components:** point `architecture`, `hmm`, `cnv`, and `evaluation` at the same versioned files as an earlier experiment and adjust the numeric parameters. No new code needed.

**New algorithm or interface:** create a new numbered file in the relevant folder (e.g. `models/hmm/02_hmm_with_ploidy.py`), then point the config at it. The old experiments remain reproducible because they still reference the original file.

## Outputs

Everything lands in the `out_dir` specified in `config.yaml`:

| File | Description |
|---|---|
| `checkpoint.pth` | Best model weights (by validation reconstruction loss) |
| `training_log.json` | Loss history per epoch |
| `latents.npy` | `(n_samples, latent_dim)` — mu vectors |
| `reconstructions.npy` | `(n_samples, n_bins)` — denormalised read-count reconstructions |
| `sample_ids.npy` | Sample ID strings matching the row order above |
| `segments.parquet` | HMM segments: `sample_id, chrom, x0, x1, cn, confidence` |
| `gene_calls.tsv` | Gene-level CN calls: `sample_id, CRT, GCH1, MDR1, PM2_PM3, crr_*` |
| `evaluation.txt` | Metrics vs Pf8 GATK ground truth (if `pf9_gt_path` set in config) |
