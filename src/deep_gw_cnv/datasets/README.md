# datasets

PyTorch `Dataset` and `DataLoader` factory code.

> **This directory contains code only.** Actual data files live under `data/`.

## Conventions

- `base.py` — abstract base class, shared split logic, common properties.
- One file per dataset variant (e.g. `gw_strain.py`, `cnv_segments.py`).
- `transforms.py` — input-space transforms shared across dataset classes.
- Register each class in `__init__.py` using the same pattern as `models/`.

## Adding a dataset

1. Subclass `BaseDataset` in a new file.
2. Document the expected directory layout in a top-of-file docstring.
3. Add to `__init__.py`:
   ```python
   from .my_dataset import MyDataset
   REGISTRY["my_dataset"] = MyDataset
   ```
4. Reference in a config: `dataset: my_dataset`
