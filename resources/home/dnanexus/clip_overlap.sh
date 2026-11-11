#!/usr/bin/env bash
# clip_overlap.sh: soft-clip the overlapping part of each read pair with fgbio ClipBam, so that a base covered by both
# reads of a pair is counted once (eggd_fgbio; DI-4247). No DNAnexus calls: the app's src/code.sh and the tests run it.
#
# usage: clip_overlap.sh -i <in.bam> -r <ref.fa.gz> -o <out_prefix> [-I <docker image>] [-t <threads>] [-m <sort mem/thread>]
#                        [-x <fgbio heap GB>] [-w <work dir>]
#   <ref.fa.gz> needs <ref.fa.gz>.fai (and .gzi if bgzipped); a sequence dictionary is made with samtools dict if absent.
# Steps (as the validation run, Confluence 4871323786 T7/T8):
#   samtools sort -n | fgbio ClipBam -c Soft --clip-overlapping-reads true | samtools sort | samtools index
# Writes <out_prefix>.clipped.bam, .clipped.bam.bai, .clipbam_metrics.txt and .provenance.tsv.
# Fails if the output does not hold exactly as many records as the input.
set -euo pipefail

IMAGE="quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0"
THREADS=$(nproc)
SORT_MEM="1G"
XMX=""
WORK=""
IN=""; REF=""; OUT=""

die() { echo "ERROR: $*" >&2; exit 1; }

while getopts "i:r:o:I:t:m:x:w:" opt; do
    case $opt in
        i) IN=$OPTARG ;;
        r) REF=$OPTARG ;;
        o) OUT=$OPTARG ;;
        I) IMAGE=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        m) SORT_MEM=$OPTARG ;;
        x) XMX=$OPTARG ;;
        w) WORK=$OPTARG ;;
        *) die "unknown option" ;;
    esac
done

if [ -z "$IN" ] || [ -z "$REF" ] || [ -z "$OUT" ]; then die "-i, -r and -o are required"; fi
[ -s "$IN" ] || die "input BAM $IN is missing or empty"
[ -s "$REF" ] || die "reference FASTA $REF is missing or empty"
[ -s "$REF.fai" ] || die "reference index $REF.fai is missing"
case $REF in *.gz) [ -s "$REF.gzi" ] || die "bgzipped reference needs $REF.gzi" ;; esac
if ! [[ $THREADS =~ ^[1-9][0-9]*$ ]]; then die "threads must be a positive integer"; fi
samtools quickcheck "$IN" || die "input $IN is not a valid BAM"
docker image inspect "$IMAGE" > /dev/null 2>&1 || die "docker image $IMAGE is not loaded"

if [ -z "$XMX" ]; then  # heap: total memory less 4 GB, between 4 and 12 GB (as the validation run)
    XMX=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 / 1024 - 4 ))
    [ "$XMX" -gt 12 ] && XMX=12
    [ "$XMX" -lt 4 ] && XMX=4
fi

OUT_DIR=$(cd "$(dirname "$OUT")" && pwd); OUT="$OUT_DIR/$(basename "$OUT")"
WORK=${WORK:-$(mktemp -d)}; mkdir -p "$WORK"; WORK=$(cd "$WORK" && pwd)
IN=$(cd "$(dirname "$IN")" && pwd)/$(basename "$IN")
REF=$(cd "$(dirname "$REF")" && pwd)/$(basename "$REF")

# fgbio (htsjdk) wants the reference, its index and a .dict side by side: link them into the work directory
case $REF in *.gz) R="$WORK/ref.fa.gz" ;; *) R="$WORK/ref.fa" ;; esac
ln -sf "$REF" "$R"; ln -sf "$REF.fai" "$R.fai"
case $REF in *.gz) ln -sf "$REF.gzi" "$R.gzi" ;; esac
dict=${REF%.gz}; dict=${dict%.*}.dict
if [ -s "$dict" ]; then ln -sf "$dict" "$WORK/ref.dict"; else samtools dict "$R" -o "$WORK/ref.dict"; fi

# the work directory and (read-only) the reference directory, at the same paths, so the links resolve in the container
REF_DIR=$(dirname "$REF")
DK=(docker run --rm -u "$(id -u):$(id -g)" -v "$WORK:$WORK" -v "$REF_DIR:$REF_DIR:ro" -w "$WORK" "$IMAGE")
t0=$(date +%s)
samtools sort -n -@ "$THREADS" -m "$SORT_MEM" -o "$WORK/qn.bam" "$IN"
"${DK[@]}" fgbio -Xmx"${XMX}"g --compression 1 ClipBam -i "$WORK/qn.bam" -o "$WORK/fg_qn.bam" -r "$R" \
    -c Soft --clip-overlapping-reads true -m "$WORK/metrics.txt"
rm -f "$WORK/qn.bam"
samtools sort -@ "$THREADS" -m "$SORT_MEM" -o "$OUT.clipped.bam" "$WORK/fg_qn.bam"
samtools index -@ "$THREADS" "$OUT.clipped.bam"
rm -f "$WORK/fg_qn.bam"
cp "$WORK/metrics.txt" "$OUT.clipbam_metrics.txt"
t1=$(date +%s)

n_in=$(samtools view -c -@ "$THREADS" "$IN")
n_out=$(samtools view -c -@ "$THREADS" "$OUT.clipped.bam")
[ "$n_in" = "$n_out" ] || die "record count changed: input $n_in, output $n_out"

{
    printf 'field\tvalue\n'
    printf 'fgbio_image\t%s\n' "$IMAGE"
    printf 'fgbio_image_id\t%s\n' "$(docker image inspect --format '{{.Id}}' "$IMAGE")"
    printf 'fgbio_version\t%s\n' "$("${DK[@]}" fgbio --version 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | awk '/Version/{print $2}')"
    printf 'samtools_version\t%s\n' "$(samtools --version | head -1 | awk '{print $2}')"
    printf 'clipbam_args\t%s\n' "-c Soft --clip-overlapping-reads true"
    printf 'threads\t%s\nsort_mem_per_thread\t%s\nfgbio_heap_gb\t%s\n' "$THREADS" "$SORT_MEM" "$XMX"
    printf 'seconds\t%s\n' "$((t1 - t0))"
    printf 'records_in\t%s\nrecords_out\t%s\n' "$n_in" "$n_out"
    samtools flagstat -@ "$THREADS" "$IN" | awk '{n=$1; $1=""; $2=""; sub(/^ +/, ""); print "flagstat_in:" $0 "\t" n}'
    samtools flagstat -@ "$THREADS" "$OUT.clipped.bam" | awk '{n=$1; $1=""; $2=""; sub(/^ +/, ""); print "flagstat_out:" $0 "\t" n}'
} > "$OUT.provenance.tsv"
echo "clip_overlap: $n_out records, $((t1 - t0)) s; wrote $OUT.clipped.bam"
