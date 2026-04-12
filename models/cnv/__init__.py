# Convenience re-export of the latest CNV caller.
# Training and wrap-up load the version specified in config.yaml via importlib.
import importlib as _il

_m = _il.import_module("cnv.01_gene_cnv_caller")

GENES_OF_INTEREST = _m.GENES_OF_INTEREST
call_gene_cnv     = _m.call_gene_cnv
call_all_genes    = _m.call_all_genes
run_cnv_calls     = _m.run_cnv_calls
