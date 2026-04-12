# Convenience re-export of the latest HMM implementation.
# Training and wrap-up load the version specified in config.yaml via importlib.
import importlib as _il

_m = _il.import_module("hmm.01_gaussian_hmm")

fit_hmm_segments    = _m.fit_hmm_segments
run_hmm_all_samples = _m.run_hmm_all_samples
