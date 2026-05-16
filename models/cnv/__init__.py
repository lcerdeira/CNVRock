# Convenience re-export of the latest CNV caller.
# Training and wrap-up load the version specified in config.yaml via importlib
# (see exp 22-35 configs: cnv: "06_gene_cnv_caller").
import importlib as _il

_m = _il.import_module("cnv.06_gene_cnv_caller")

# Re-export only the symbols that exist on the current KpSC caller.
# Older Pf-era helpers (GENES_OF_INTEREST, call_gene_cnv, call_all_genes) are
# no longer part of the active API; they have been removed.
run_cnv_calls = _m.run_cnv_calls
