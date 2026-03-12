module load HGI/common/nextflow/25.04.6

bsub \
    -J read-counts \
    -o log.o \
    -e log.e \
    -n 1 \
    -G gsu-pipelines \
    -M 8000 \
    -R 'select[mem>=8000] rusage[mem=8000] span[hosts=1]' \
    -q basement \
    "nextflow run main.nf -profile sanger_lsf -resume"
