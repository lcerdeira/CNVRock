# deep-gw-cnv

Convolutional VAE for genome-wide CNV detection from read-count data.

```
TODO:


```

## Layout

```
data/         zarr stores (read counts) and results
assets/       sample manifests and BAM/CRAM path lists
setup/        scripts to extract read counts and convert to zarr
models/       all model code
  train.py          entry point — run any experiment
  architectures/    model definitions (ConvVAE)
  training/         dataset, trainer, inference helpers
  experiments/      one folder per experiment, each self-contained
```

## Running an experiment

```bash
cd models/experiments/01
bash run_mac.sh        # local (MPS or CPU)
bash run_cluster.sh    # LSF GPU cluster
```

Or directly:
```bash
python3 models/train.py models/experiments/01/config.yaml
```

## Adding a new experiment

```bash
cp -r models/experiments/01 models/experiments/02
# edit models/experiments/02/config.yaml
```

Outputs (`latents.npy`, `reconstructions.npy`, `sample_ids.npy`, `checkpoint.pth`) are written to the `out_dir` defined in the config.

## Setup

See [setup/](setup/) for Nextflow pipelines that extract read counts from BAMs/CRAMs and write zarr stores.

## Quick start
