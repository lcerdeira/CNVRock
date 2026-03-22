module load HGI/common/nextflow/25.04.6

bsub \
    -J "rc2zarr" \
    -o "log.o" \
    -e "log.e" \
    -q basement \
    -G team342 \
    -R "select[mem>48000] rusage[mem=48000]" \
    -M 48000 \
    -n 8 \
    "PYTHONUNBUFFERED=1 python3 -u readcounts_to_zarr.py"