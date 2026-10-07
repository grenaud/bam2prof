#!/bin/bash
# Checks that every way bam2prof can load the -fa reference gives the same profiles, and reports peak memory:
#   sorted+indexed BAM          -> reference mapped 10Mb at a time (windowed, with prefetch)
#   unsorted BAM, no index      -> read sequentially, whole reference mapped
#   same, piped through stdin   -> likewise
#   sorted BAM piped via stdin  -> sequential read, windowed reference
# The sorted+indexed run is the reference answer; the other three must match it exactly.
#
# Usage: ./test_reference_modes.sh [-b UNSORTED_BAM] [-n NREADS] [-r REF] [-o OUTDIR]
#   -b  unsorted (e.g. queryname-sorted) BAM, default testData/chr2_name.bam
#   -n  use only its first NREADS records (default 2000000; 0 = all, slow for the 5.5GB file)
#   -r  reference fasta (default testData/reference/chr2_only.fa)
#   -o  output dir (default results/refmodes_<timestamp>)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
BAM=testData/chr2_name.bam; N=2000000; REF=testData/reference/chr2_only.fa; OUT=""
while getopts "b:n:r:o:h" o; do case $o in b) BAM=$OPTARG;; n) N=$OPTARG;; r) REF=$OPTARG;; o) OUT=$OPTARG;; *) sed -n '2,13p' "$0"; exit 1;; esac; done
[ -n "$OUT" ] || OUT=results/refmodes_$(date +%Y%m%d_%H%M%S)
SAM=lib/samtools/samtools; B2P=./src/bam2prof; ARGS="-classic -comp -around 12 -length 20 -fa $REF"
mkdir -p "$OUT"

if [ "$N" -gt 0 ]; then
    { $SAM view -h "$BAM" || true; } | head -n $((N + $($SAM view -H "$BAM" | wc -l))) | $SAM view -b -o "$OUT/unsorted.bam" -
else
    ln -sf "$(realpath "$BAM")" "$OUT/unsorted.bam"
fi
$SAM sort -o "$OUT/sorted.bam" "$OUT/unsorted.bam"; $SAM index "$OUT/sorted.bam"

run(){ # name, command reading from $1 (a file or "-")
    local name=$1; shift
    /usr/bin/time -v "$@" 2> "$OUT/$name.err" >/dev/null
}
run sorted_indexed  $B2P $ARGS -o "$OUT/sorted_indexed"  "$OUT/sorted.bam"
run unsorted_file   $B2P $ARGS -o "$OUT/unsorted_file"   "$OUT/unsorted.bam"
cat "$OUT/unsorted.bam"           | run unsorted_pipe  $B2P $ARGS -o "$OUT/unsorted_pipe" -
$SAM view -b "$OUT/sorted.bam"    | run sorted_pipe    $B2P $ARGS -o "$OUT/sorted_pipe" -

# profile files named after the input; compare contents only
norm(){ for f in "$OUT/$1"/*; do echo "$(basename "$f" | sed -E 's/^.*_classic/X_classic/') $(md5sum < "$f")"; done | sort; }
fail=0
printf "%-16s %-10s %s\n" run "peak RSS" "vs sorted_indexed"
for d in sorted_indexed unsorted_file unsorted_pipe sorted_pipe; do
    rss=$(grep 'Maximum resident' "$OUT/$d.err" | awk '{printf "%.0f MB", $NF/1024}')
    if diff <(norm "$d") <(norm sorted_indexed) >/dev/null; then res=identical; else res=DIFFERENT; fail=1; fi
    printf "%-16s %-10s %s\n" "$d" "$rss" "$res"
done
exit $fail
