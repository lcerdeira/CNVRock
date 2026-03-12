# notebooks

Jupyter notebooks for exploratory data analysis, post-hoc experiment
analysis, and interactive demos. Notebooks are **not** part of the training
pipeline — all reusable logic should be extracted into `src/deep_gw_cnv/`.

## Naming convention

```
NNN_short_description.ipynb
```

e.g. `001_data_exploration.ipynb`, `002_loss_curve_comparison.ipynb`

## Conventions

- Clear all cell outputs before committing (`Kernel → Restart & Clear Output`)
  to keep diffs readable.
- If a notebook develops into a reusable utility, move the logic to
  `src/deep_gw_cnv/utils/` and call it from the notebook.
- Notebooks should be self-contained: include data-loading cells at the top
  so a fresh kernel can run the whole thing top-to-bottom.
