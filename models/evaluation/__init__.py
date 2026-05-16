# Convenience re-export of the latest evaluation module.
# Training and wrap-up load the version specified in config.yaml via importlib.
import importlib as _il

_m = _il.import_module("evaluation.04_kpsc_evaluation")

run_evaluation = _m.run_evaluation
