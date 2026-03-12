# augmentation

Signal-space transforms applied during training. These are library classes,
not scripts — they are imported by `Dataset` subclasses in `datasets/`.

## Conventions

- Each transform exposes a `__call__(self, sample)` interface compatible with
  `torchvision.transforms.Compose`.
- Transforms are **stateless** with respect to the training loop; all
  randomness comes from seeded RNGs passed at construction or drawn per call.
- Parameters come from the `augmentation` block in the experiment config; the
  transform classes should accept a `cfg` dict and ignore unknown keys.

## Adding a transform

1. Create or add to an existing file in this directory.
2. Ensure the class is exported from `__init__.py`.
3. Reference it in the `augmentation` block of a config.
