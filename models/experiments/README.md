# experiments/

Each sub-folder is a self-contained experiment.

```
01/
  config.yaml       all hyperparameters and paths
  run_mac.sh        run locally (MPS or CPU)
  run_cluster.sh    submit to LSF GPU cluster
  experiment.py     thin wrapper (calls train.py) — prefer the shell scripts
```

**To add experiment 02:**
```bash
cp -r 01 02
# edit 02/config.yaml — at minimum change out_dir
```

Outputs land in the `out_dir` specified in `config.yaml`:
- `checkpoint.pth` — best model weights (by reconstruction loss)
- `latents.npy` — `(n_samples, latent_dim)` mu vectors
- `reconstructions.npy` — `(n_samples, n_bins)` in log2(count+1) space
- `sample_ids.npy` — matching sample ID strings
