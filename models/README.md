# models/

## Entry point

```bash
python3 train.py path/to/experiment/config.yaml
```

Loads config, builds the model, trains, then runs inference on every sample.

## Modules

| Path | What it does |
|---|---|
| `train.py` | Reads a `config.yaml` and runs train + inference |
| `architectures/01_conv_vae.py` | `ConvVAE` — 1-D convolutional encoder/decoder |
| `training/dataset.py` | `ReadCountDataset` — reads a zarr store, log2-normalises |
| `training/trainer.py` | `train_vae` — training loop with KL annealing + early stopping |
| `training/wrap_up.py` | `run_inference` + `run_hmm_all_samples` — inference, HMM segmentation, saves all outputs |
| `experiments/` | One folder per experiment, each with a `config.yaml` and shell scripts |
