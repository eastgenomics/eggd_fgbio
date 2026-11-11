"""Tests for src/code.sh, the DNAnexus entry point, with DNAnexus and Docker commands replaced by stubs.

The stubs stand in for dx-download-all-inputs, dx upload, dx-jobutil-add-output, docker and ~/clip_overlap.sh (the
pipeline, tested in test_clip_overlap.py), so these tests cover what the entry point itself does: loading the image
input, choosing the tagged image, passing the image and version to the pipeline, and recording the image file in the
provenance. HOME is a temporary directory, as /home/dnanexus is on the worker.
"""
import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

CODE = Path(__file__).resolve().parents[1] / "src/code.sh"
IMAGE = "quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0"

STUB_DOCKER = """#!/bin/bash
# docker stub: `docker load` prints $DOCKER_LOAD_OUT; nothing else is used by the entry point
[ "$1" = load ] && printf '%b\\n' "$DOCKER_LOAD_OUT"
exit 0
"""
STUB_DX = """#!/bin/bash
# dx stub: `dx upload <path> --wait --brief` prints a file ID made from the file name
[ "$1" = upload ] && echo "file-$(basename "$2")"
exit 0
"""
STUB_ADD_OUTPUT = """#!/bin/bash
echo "$*" >> "$HOME/added_outputs.txt"
"""
STUB_PIPELINE = """#!/bin/bash
# clip_overlap.sh stub: records its arguments and writes the four outputs for -o <prefix>
echo "$*" > "$HOME/pipeline_args.txt"
while getopts "i:r:o:I:V:w:" opt; do [ "$opt" = o ] && out=$OPTARG; done
for s in clipped.bam clipped.bam.bai clipbam_metrics.txt; do echo x > "$out.$s"; done
printf 'field\\tvalue\\tqc_failed\\nfgbio_version\\t4.1.1\\t\\n' > "$out.provenance.tsv"
"""


def write_exe(path, text):
    """Write an executable stub script."""
    path.write_text(text)
    path.chmod(0o755)


@pytest.fixture
def worker(tmp_path):
    """Make a fake worker: HOME with inputs, the pipeline stub and executable description, and stubs on PATH.

    Returns (environment, path of the image file). The input variables are those dx-download-all-inputs would set.
    """
    home = tmp_path / "home"
    (home / "in").mkdir(parents=True)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_exe(bin_dir / "docker", STUB_DOCKER)
    write_exe(bin_dir / "dx", STUB_DX)
    write_exe(bin_dir / "dx-jobutil-add-output", STUB_ADD_OUTPUT)
    write_exe(bin_dir / "dx-download-all-inputs", "#!/bin/bash\nexit 0\n")
    write_exe(bin_dir / "samtools", "#!/bin/bash\necho 'samtools 1.24'\n")
    write_exe(home / "clip_overlap.sh", STUB_PIPELINE)
    (home / "dnanexus-executable.json").write_text(json.dumps({"name": "eggd_fgbio", "version": "1.0.0", "id": "app-TEST"}))
    files = {"bam": "s.bam", "ref": "ref.fa.gz", "fai": "ref.fa.gz.fai", "gzi": "ref.fa.gz.gzi",
             "docker": "fgbio.docker.tar.gz"}
    for name in files.values():
        (home / "in" / name).write_bytes(name.encode())
    env = dict(os.environ, HOME=str(home), PATH=f"{bin_dir}:{os.environ['PATH']}",
               bam_path=str(home / "in/s.bam"), bam_prefix="s",
               reference_fasta_path=str(home / "in/ref.fa.gz"), reference_fasta_name="ref.fa.gz",
               reference_fai_path=str(home / "in/ref.fa.gz.fai"), reference_fai_name="ref.fa.gz.fai",
               reference_gzi='{"$dnanexus_link": "file-GZI"}', reference_gzi_path=str(home / "in/ref.fa.gz.gzi"),
               reference_gzi_name="ref.fa.gz.gzi",
               fgbio_docker='{"$dnanexus_link": {"project": "project-P", "id": "file-DOCKER"}}',
               fgbio_docker_path=str(home / "in/fgbio.docker.tar.gz"), fgbio_docker_name="fgbio.docker.tar.gz",
               fgbio_version="4.1.1", output_prefix="", DOCKER_LOAD_OUT=f"Loaded image: {IMAGE}")
    return env, home / "in/fgbio.docker.tar.gz"


def run_main(env):
    """Source code.sh and run main() as DNAnexus would, returning the completed process."""
    return subprocess.run(["bash", "-c", f'source "{CODE}"; main'], env=env, capture_output=True, text=True)


def test_one_tagged_image_is_run_with_the_given_version(worker):
    """Check that the loaded image's tag and the fgbio_version input are passed to the pipeline."""
    env, _ = worker
    r = run_main(env)
    assert r.returncode == 0, r.stderr
    args = (Path(env["HOME"]) / "pipeline_args.txt").read_text().split()
    assert args[args.index("-I") + 1] == IMAGE and args[args.index("-V") + 1] == "4.1.1"


def test_provenance_records_the_image_file_and_app(worker):
    """Check that the provenance gains the image file name and ID, its SHA-256 and the app, and outputs are declared."""
    env, docker_file = worker
    sha = hashlib.sha256(docker_file.read_bytes()).hexdigest()
    r = run_main(env)
    assert r.returncode == 0, r.stderr
    rows = [line.split("\t") for line in (Path(env["HOME"]) / "out/s.provenance.tsv").read_text().splitlines()]
    prov = {row[0]: row[1] for row in rows[1:]}
    assert prov["fgbio_docker"] == "fgbio.docker.tar.gz file-DOCKER"
    assert prov["fgbio_docker_sha256"] == sha
    assert prov["app"] == "eggd_fgbio 1.0.0 app-TEST"
    added = (Path(env["HOME"]) / "added_outputs.txt").read_text().splitlines()
    assert [line.split()[0] for line in added] == ["clipped_bam", "clipped_bai", "clipbam_metrics", "provenance"]


@pytest.mark.parametrize("load_out", [
    "Loaded image ID: sha256:c2b12e4102d185d52c0a7f8dae3725a4cbf4cdc3e4ec091bb959c29eeaa15862",  # untagged image
    f"Loaded image: {IMAGE}\\nLoaded image: quay.io/biocontainers/samtools:1.24",              # two tagged images
    "",                                                                                       # nothing loaded
])
def test_image_file_without_exactly_one_tag_fails(worker, load_out):
    """Check that an image file with no tagged image, or more than one, fails before the pipeline or any upload."""
    env, _ = worker
    env["DOCKER_LOAD_OUT"] = load_out
    r = run_main(env)
    assert r.returncode == 1 and "ERROR: fgbio_docker (fgbio.docker.tar.gz) must hold exactly one tagged image" in r.stderr
    assert not (Path(env["HOME"]) / "pipeline_args.txt").exists()
    assert not (Path(env["HOME"]) / "added_outputs.txt").exists()
