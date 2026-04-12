# Convenience re-export of the latest architecture.
# Training and wrap-up load the version specified in config.yaml via importlib.
import importlib as _il

_m = _il.import_module("architectures.01_conv_vae")

ConvVAE       = _m.ConvVAE
N_BINS_RAW    = _m.N_BINS_RAW
N_BINS_PADDED = _m.N_BINS_PADDED
