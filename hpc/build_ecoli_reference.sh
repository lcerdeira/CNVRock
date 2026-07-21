#!/bin/bash
#SBATCH --job-name=ecoli_ref
#SBATCH --output=/home/lshlt19/CNVRock/logs/ecoli_ref_%j.out
#SBATCH --error=/home/lshlt19/CNVRock/logs/ecoli_ref_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Index the S. aureus extended reference and emit the 1 kb interval list.
#
# The interval list is generated FROM the sequence dictionary rather than
# from the FASTA directly, so the @SQ records (name, length, M5) are the ones
# GATK will itself compute. Hand-rolling the header risks an M5 mismatch that
# only surfaces once CollectReadCounts runs, i.e. after the whole array has
# already spent its download time.
#
# Usage:
#   sbatch hpc/build_ecoli_reference.sh

set -euo pipefail

REPO_DIR="/home/lshlt19/CNVRock"
ENV_BIN="/home/lshlt19/miniconda3/envs/cnvrock/bin"
export PATH="$ENV_BIN:/home/lshlt19/miniconda3/bin:$PATH"

module load bwa/0.718 samtools/1.20 gatk/4.6.0.0 java/20.0.1

REF="$REPO_DIR/assets/ecoli_ref/EC958_extended.fasta"
DICT="${REF%.fasta}.dict"
INTERVALS="$REPO_DIR/assets/ecoli_ref/EC958_extended_1kb.interval_list"

[[ -f "$REF" ]] || { echo "missing reference: $REF" >&2; exit 1; }

echo "== bwa index =="
bwa index "$REF"

echo "== samtools faidx =="
samtools faidx "$REF"

echo "== sequence dictionary =="
rm -f "$DICT"
gatk CreateSequenceDictionary -R "$REF" -O "$DICT"

echo "== 1 kb interval list =="
python3 - "$DICT" "$INTERVALS" <<'PY'
import sys
dict_path, out_path = sys.argv[1], sys.argv[2]
BIN = 1000
header, contigs = [], []
with open(dict_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("@HD"):
            header.append("@HD\tVN:1.6\tSO:unsorted")
        elif line.startswith("@SQ"):
            header.append(line)
            f = dict(p.split(":", 1) for p in line.split("\t")[1:] if ":" in p)
            contigs.append((f["SN"], int(f["LN"])))
if not any(h.startswith("@HD") for h in header):
    header.insert(0, "@HD\tVN:1.6\tSO:unsorted")

n = 0
with open(out_path, "w") as out:
    for h in header:
        out.write(h + "\n")
    for name, length in contigs:
        for start in range(1, length + 1, BIN):
            end = min(start + BIN - 1, length)
            out.write(f"{name}\t{start}\t{end}\t+\t.\n")
            n += 1
print(f"  {len(contigs)} contigs, {n} intervals -> {out_path}")
PY

echo "== done =="
ls -la "$REF".* "$DICT" "$INTERVALS"
