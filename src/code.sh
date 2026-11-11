#!/bin/bash
# Input variables (bam_path, reference_fasta_name, ...) are set by DNAnexus and dx-download-all-inputs.
# shellcheck disable=SC2154
# eggd_fgbio: fgbio tools on DNAnexus. v1.0.0 runs fgbio 4.1.1 ClipBam to soft-clip overlapping read pairs in a BAM,
# the input to AMBER in the Atlas CNV workflow's PURPLE arm (DI-4247). The pipeline itself is
# /home/dnanexus/clip_overlap.sh; this file only handles DNAnexus input and output.
set -exo pipefail

# fgbio image, bundled with the app (runSpec.bundledDepends) as a `docker save` tar of
# quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace
IMAGE="quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0"
IMAGE_TAR="/home/dnanexus/fgbio_image/fgbio_4.1.1--hdfd78af_0.docker.tar"
IMAGE_TAR_SHA256="6740248d07c7049def1331f030f27e35cc4096c46cce1913021dbe83f6bb99c6"

main() {
    echo "=== Checking and loading the fgbio image"
    echo "${IMAGE_TAR_SHA256}  ${IMAGE_TAR}" | sha256sum -c -
    docker load -i "$IMAGE_TAR"
    docker image inspect --format '{{.Id}}' "$IMAGE"
    rm -f "$IMAGE_TAR"

    echo "=== Downloading inputs"
    dx-download-all-inputs --parallel

    # The reference FASTA and its indexes must be named as one set: <fasta>, <fasta>.fai and (bgzipped) <fasta>.gzi.
    mkdir -p ~/ref
    mv "$reference_fasta_path" ~/ref/
    ref=~/ref/"$reference_fasta_name"
    [ "$reference_fai_name" = "${reference_fasta_name}.fai" ] || \
        { echo "ERROR: reference_fai ($reference_fai_name) is not ${reference_fasta_name}.fai" >&2; exit 1; }
    mv "$reference_fai_path" ~/ref/
    if [ -n "$reference_gzi" ]; then
        [ "$reference_gzi_name" = "${reference_fasta_name}.gzi" ] || \
            { echo "ERROR: reference_gzi ($reference_gzi_name) is not ${reference_fasta_name}.gzi" >&2; exit 1; }
        mv "$reference_gzi_path" ~/ref/
    fi

    prefix="${output_prefix:-${bam_prefix}}"
    case $prefix in
        ""|*/*|.*) echo "ERROR: output prefix '$prefix' must be a plain file name" >&2; exit 1 ;;
    esac

    echo "=== Clipping overlapping read pairs"
    samtools --version | head -1
    mkdir -p ~/out
    bash /home/dnanexus/clip_overlap.sh -i "$bam_path" -r "$ref" -o ~/out/"$prefix" -I "$IMAGE" -w ~/work
    printf 'app\t%s\n' "$(python3 -c 'import json; d = json.load(open("/home/dnanexus/dnanexus-executable.json")); print(d["name"], d.get("version", ""), d["id"])' 2>/dev/null || echo unknown)" \
        >> ~/out/"$prefix".provenance.tsv

    echo "=== Uploading outputs"
    for pair in "clipped_bam:clipped.bam" "clipped_bai:clipped.bam.bai" "clipbam_metrics:clipbam_metrics.txt" \
                "provenance:provenance.tsv"; do
        id=$(dx upload ~/out/"$prefix.${pair#*:}" --wait --brief)
        dx-jobutil-add-output "${pair%%:*}" "$id" --class=file
    done
}
