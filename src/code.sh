#!/bin/bash
# shellcheck disable=SC2154  # input variables (bam_path, reference_fasta_name, ...) are set by DNAnexus and dx-download-all-inputs
# eggd_fgbio: fgbio tools on DNAnexus. v1.0.0 runs fgbio ClipBam to soft-clip overlapping read pairs in a BAM, the input
# to AMBER in the Atlas CNV workflow's PURPLE arm (DI-4247). The fgbio Docker image and the fgbio version it must hold are
# inputs, so a new image can be used without changing the app; the job fails if the two do not match. The pipeline itself
# is ~/clip_overlap.sh (HOME is /home/dnanexus on the worker); this file only handles DNAnexus input and output.
set -exo pipefail

main() {
    echo "=== Downloading inputs"
    dx-download-all-inputs --parallel

    echo "=== Loading the fgbio image"
    image_sha256=$(sha256sum "$fgbio_docker_path" | cut -d' ' -f1)
    load_out=$(docker load -i "$fgbio_docker_path")
    echo "$load_out"
    # exactly one tagged image: the tag is what the pipeline runs
    image=$(sed -n 's/^Loaded image: //p' <<< "$load_out")
    if [ -z "$image" ] || [ "$(wc -l <<< "$image")" -ne 1 ]; then
        echo "ERROR: fgbio_docker ($fgbio_docker_name) must hold exactly one tagged image" >&2; exit 1
    fi
    rm -f "$fgbio_docker_path"
    # shellcheck disable=SC2016  # "$dnanexus_link" is a JSON key for Python, not a shell expansion
    docker_id=$(python3 -c 'import json, sys; v = json.loads(sys.argv[1])["$dnanexus_link"]; print(v if isinstance(v, str) else v["id"])' "$fgbio_docker")

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
    bash ~/clip_overlap.sh -i "$bam_path" -r "$ref" -o ~/out/"$prefix" -I "$image" -V "$fgbio_version" -w ~/work
    # the executable that ran (name, version, ID) is part of the provenance; failing to read any of it fails the job.
    # Apps have a version; development applets do not, and are recorded as "applet".
    app_desc=$(python3 -c '
import json, os, sys
d = json.load(open(os.path.expanduser("~/dnanexus-executable.json")))
name, eid = d.get("name", ""), d.get("id", "")
ver = d.get("version", "") if eid.startswith("app-") else "applet"
if not (name and eid and ver):
    sys.exit("ERROR: the executable description lacks a name, ID or app version")
print(name, ver, eid)')
    {
        printf 'app\t%s\t\n' "$app_desc"
        printf 'fgbio_docker\t%s %s\t\n' "$fgbio_docker_name" "$docker_id"
        printf 'fgbio_docker_sha256\t%s\t\n' "$image_sha256"
    } >> ~/out/"$prefix".provenance.tsv

    echo "=== Uploading outputs"
    for pair in "clipped_bam:clipped.bam" "clipped_bai:clipped.bam.bai" "clipbam_metrics:clipbam_metrics.txt" \
                "provenance:provenance.tsv"; do
        id=$(dx upload ~/out/"$prefix.${pair#*:}" --wait --brief)
        dx-jobutil-add-output "${pair%%:*}" "$id" --class=file
    done
}
