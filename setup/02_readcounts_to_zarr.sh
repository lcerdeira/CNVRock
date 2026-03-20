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
    "python3 readcounts_to_zarr.py"