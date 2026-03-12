# data

All data files. **Fully gitignored** — nothing in this directory is tracked
by version control.

## Layout

```
data/
├── raw/          ← original, immutable source files; never modified after download
└── processed/    ← outputs of scripts/preprocess.py; can always be regenerated
```

## Acquiring raw data

[Describe where to download the data or how to request access, e.g.:]
```bash
# Example — replace with the actual command or URL
wget -P data/raw/ https://example.org/dataset.tar.gz
tar -xf data/raw/dataset.tar.gz -C data/raw/
```

## Reproducing processed data

Once raw data is in place, run the preprocessing pipeline once:
```bash
python scripts/preprocess.py --raw data/raw/ --out data/processed/
```

Expected disk usage: raw ≈ ? GB, processed ≈ ? GB.

## .gitkeep files

The `raw/` and `processed/` subdirectories contain `.gitkeep` placeholder
files so the directory structure is tracked even when empty.
