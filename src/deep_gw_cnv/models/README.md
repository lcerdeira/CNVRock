# models

Architecture definitions. Each model is a `torch.nn.Module` in its own file
and registered in `__init__.py` so configs can reference models by string name.

## Adding a model

1. Create `src/deep_gw_cnv/models/my_model.py`
2. Define a class that inherits from `nn.Module`
3. Register it in `__init__.py`:
   ```python
   from .my_model import MyModel
   REGISTRY["my_model"] = MyModel
   ```
4. Reference it in a config: `model: my_model`

## Conventions

- Constructors accept a single `cfg` dict (no bare positional hyperparams).
- Keep `forward()` free of side effects — this makes tracing and scripting
  straightforward.
- Prefer composition over inheritance; keep blocks in separate files if they
  are reused across architectures.
