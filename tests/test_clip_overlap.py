"""Tests for resources/home/dnanexus/clip_overlap.sh on a small synthetic BAM.

Needs samtools and bgzip on PATH and the fgbio image (quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0) in the local docker
daemon; CI builds samtools/htslib 1.24 (the version the app uses) and pulls the image by digest. The reference and reads
are generated here (random sequence, fixed seed), so no patient data.
"""
import os
import random
import re
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "resources/home/dnanexus/clip_overlap.sh"
IMAGE = "quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0"
FGBIO_VERSION = "4.1.1"
SCRIPT_CMD = ["bash", str(SCRIPT), "-I", IMAGE, "-V", FGBIO_VERSION]  # the script with the image and its fgbio version
READ_LEN = 100
# name: (pos1, pos2), 1-based leftmost positions of read 1 (forward) and read 2 (reverse) on chr1
PAIRS = {
    "overlap_40": (1001, 1061),       # reads overlap by 40 bases
    "overlap_full": (2001, 2001),     # reads cover the same 100 bases
    "abutting": (3001, 3101),         # no overlap: read 2 starts after read 1 ends
    "apart": (4001, 4301),            # far apart
}


def run(cmd, **kw):
    """Run a command, failing the test on a non-zero exit, and return the completed process."""
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def revcomp(s):
    """Return the reverse complement of a DNA sequence."""
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    """Make the synthetic reference and BAM used by every test.

    The reference is two random 6 kb contigs (fixed seed), written plain and bgzipped with their indexes. The BAM holds
    the read pairs in PAIRS on chr1, all FR orientation with 100-base reads, and a pair with an unmapped mate on chr2.
    Returns the directory holding ref.fa and ref.fa.gz with their indexes, and in.bam (coordinate-sorted, not indexed).
    """
    d = tmp_path_factory.mktemp("data")
    rnd = random.Random(4247)
    seq = {c: "".join(rnd.choice("ACGT") for _ in range(6000)) for c in ("chr1", "chr2")}
    fa = d / "ref.fa"
    fa.write_text("".join(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n" for c, s in seq.items()))
    run(["samtools", "faidx", str(fa)])
    with open(d / "ref.fa.gz", "wb") as gz:  # bgzipped copy, as the production reference
        subprocess.run(["bgzip", "-c", str(fa)], check=True, stdout=gz)
    run(["samtools", "faidx", str(d / "ref.fa.gz")])
    lines = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{c}\tLN:{len(s)}" for c, s in seq.items()] + \
            ["@RG\tID:rg1\tSM:synthetic\tPL:ILLUMINA"]
    q = "I" * READ_LEN
    for name, (p1, p2) in PAIRS.items():
        tlen = p2 + READ_LEN - p1
        r1 = seq["chr1"][p1 - 1:p1 - 1 + READ_LEN]
        r2 = seq["chr1"][p2 - 1:p2 - 1 + READ_LEN]
        lines.append(f"{name}\t99\tchr1\t{p1}\t60\t{READ_LEN}M\t=\t{p2}\t{tlen}\t{r1}\t{q}\tRG:Z:rg1")
        lines.append(f"{name}\t147\tchr1\t{p2}\t60\t{READ_LEN}M\t=\t{p1}\t{-tlen}\t{r2}\t{q}\tRG:Z:rg1")
    r = seq["chr2"][499:599]
    lines.append(f"mate_unmapped\t73\tchr2\t500\t60\t{READ_LEN}M\t=\t500\t0\t{r}\t{q}\tRG:Z:rg1")
    lines.append(f"mate_unmapped\t133\tchr2\t500\t0\t*\t=\t500\t0\t{revcomp(r)}\t{q}\tRG:Z:rg1")
    sam = d / "in.sam"
    sam.write_text("\n".join(lines) + "\n")
    run(["samtools", "sort", "-o", str(d / "in.bam"), str(sam)])
    return d


def records(bam):
    """Return the records of a BAM keyed by read name and read-1/read-2 bits.

    Each value holds flag, pos, cigar, seq and qual. Tags are left out because fgbio may update them.
    """
    out = {}
    for line in run(["samtools", "view", str(bam)]).stdout.splitlines():
        f = line.split("\t")
        out[(f[0], int(f[1]) & 0xC0)] = {"flag": f[1], "pos": int(f[3]), "cigar": f[5], "seq": f[9], "qual": f[10]}
    return out


def ref_end(pos, cigar):
    """Return the 1-based reference position of the last base an alignment covers."""
    return pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X") - 1


@pytest.fixture(scope="module")
def clipped(data, tmp_path_factory):
    """Run clip_overlap.sh on the synthetic BAM with the plain reference and return its output directory."""
    o = tmp_path_factory.mktemp("out")
    run([*SCRIPT_CMD, "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"), "-o", str(o / "s"), "-t", "2",
         "-m", "100M", "-x", "1", "-w", str(o / "work")])
    return o


@pytest.fixture(scope="module")
def clipped_bgzf(data, tmp_path_factory):
    """Run clip_overlap.sh with the bgzipped reference (.fa.gz, .fai, .gzi), as in production; return its output directory."""
    o = tmp_path_factory.mktemp("out_bgzf")
    run([*SCRIPT_CMD, "-i", str(data / "in.bam"), "-r", str(data / "ref.fa.gz"), "-o", str(o / "s"), "-t", "2",
         "-m", "100M", "-x", "1", "-w", str(o / "work")])
    return o


def test_outputs_exist(clipped):
    """Check that the clipped BAM, its index, the ClipBam metrics and the provenance are written and not empty."""
    for suffix in (".clipped.bam", ".clipped.bam.bai", ".clipbam_metrics.txt", ".provenance.tsv"):
        assert (clipped / f"s{suffix}").stat().st_size > 0, suffix


def test_record_count_unchanged(data, clipped):
    """Check that clipping neither adds nor removes records, counted with samtools view -c."""
    n_in, n_out = (int(run(["samtools", "view", "-c", str(bam)]).stdout) for bam in (data / "in.bam", clipped / "s.clipped.bam"))
    assert n_out == n_in == 2 * len(PAIRS) + 2


def test_overlapping_pairs_no_longer_overlap(clipped):
    """Check that both overlapping pairs are soft-clipped and no longer overlap, and that the 40-base overlap is split.

    fgbio clips each mate back to the midpoint of the overlap, so each read of the 40-base pair keeps 80 aligned bases.
    """
    rec = records(clipped / "s.clipped.bam")
    for name in ("overlap_40", "overlap_full"):
        r1, r2 = rec[(name, 0x40)], rec[(name, 0x80)]
        p1, c1, p2, c2 = r1["pos"], r1["cigar"], r2["pos"], r2["cigar"]
        assert "S" in c1 + c2, (name, c1, c2)
        assert ref_end(p1, c1) < p2 or ref_end(p2, c2) < p1, (name, p1, c1, p2, c2)
    # the 40 overlapping bases are split evenly between the two reads: each keeps 100 - 20 aligned bases
    r1, r2 = rec[("overlap_40", 0x40)], rec[("overlap_40", 0x80)]
    assert [ref_end(r["pos"], r["cigar"]) - r["pos"] + 1 for r in (r1, r2)] == [READ_LEN - 20, READ_LEN - 20], (r1, r2)


def test_clipping_keeps_flags_sequence_and_qualities(data, clipped):
    """Check that soft clipping keeps every base: flag, sequence and qualities are unchanged for all reads."""
    before, after = records(data / "in.bam"), records(clipped / "s.clipped.bam")
    for key in before:
        for field in ("flag", "seq", "qual"):
            assert after[key][field] == before[key][field], (key, field)


def test_other_reads_keep_position_and_cigar(data, clipped):
    """Check that reads not overlapping their mate keep flag, position, CIGAR, sequence and qualities (tags may change)."""
    before, after = records(data / "in.bam"), records(clipped / "s.clipped.bam")
    for key in before:
        if key[0] not in ("overlap_40", "overlap_full"):
            assert after[key] == before[key], key


def test_bgzipped_reference_gives_the_same_records(clipped, clipped_bgzf):
    """Check that the bgzipped and plain references give the same flags, positions, CIGARs, sequences and qualities."""
    assert records(clipped_bgzf / "s.clipped.bam") == records(clipped / "s.clipped.bam")


def provenance(path):
    """Read a provenance TSV, check it has three columns throughout, and return field -> (value, qc_failed)."""
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    assert rows[0] == ["field", "value", "qc_failed"] and all(len(r) == 3 for r in rows), rows
    return {r[0]: (r[1], r[2]) for r in rows[1:]}


def test_provenance(clipped):
    """Check the provenance fields, versions, image ID, record counts and flagstat rows of a run."""
    prov = provenance(clipped / "s.provenance.tsv")
    assert prov["fgbio_version"][0] == prov["fgbio_version_expected"][0] == FGBIO_VERSION
    assert prov["fgbio_image"][0] == IMAGE
    assert prov["fgbio_image_id"][0].startswith("sha256:")
    assert re.fullmatch(r"\d+\.\d+(\.\d+)?", prov["samtools_version"][0])
    assert prov["records_in"] == prov["records_out"] == (str(2 * len(PAIRS) + 2), "")
    # flagstat keys are fixed descriptions with separate QC-passed and QC-failed counts
    for side in ("in", "out"):
        assert prov[f"flagstat_{side}:total (QC-passed reads + QC-failed reads)"] == (str(2 * len(PAIRS) + 2), "0")
        assert prov[f"flagstat_{side}:paired in sequencing"] == (str(2 * len(PAIRS) + 2), "0")
        assert not any("%" in k for k in prov if k.startswith(f"flagstat_{side}:"))


def test_unreadable_provenance_value_fails(data, tmp_path):
    """Check that a samtools whose --version prints nothing fails the run instead of writing incomplete provenance."""
    real = shutil.which("samtools")
    fake = tmp_path / "bin" / "samtools"
    fake.parent.mkdir()
    fake.write_text(f'#!/bin/bash\n[ "$1" = --version ] && exit 0\nexec {real} "$@"\n')
    fake.chmod(0o755)
    env = dict(os.environ, PATH=f"{fake.parent}:{os.environ['PATH']}")
    r = subprocess.run([*SCRIPT_CMD, "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"), "-o",
                        str(tmp_path / "s"), "-t", "2", "-m", "100M", "-x", "1", "-w", str(tmp_path / "work")],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 1 and "provenance value samtools_version could not be read" in r.stderr
    assert not (tmp_path / "s.provenance.tsv").exists()


def test_missing_fai_fails(data, tmp_path):
    """Check that a reference without its .fai index fails with an ERROR: message."""
    (tmp_path / "ref.fa").write_text((data / "ref.fa").read_text())
    r = subprocess.run([*SCRIPT_CMD, "-i", str(data / "in.bam"), "-r", str(tmp_path / "ref.fa"),
                        "-o", str(tmp_path / "s")], capture_output=True, text=True)
    assert r.returncode == 1 and "ERROR: reference index" in r.stderr


def test_not_a_bam_fails(data, tmp_path):
    """Check that an input that is not a BAM fails with an ERROR: message."""
    bad = tmp_path / "bad.bam"
    bad.write_text("not a bam\n")
    r = subprocess.run([*SCRIPT_CMD, "-i", str(bad), "-r", str(data / "ref.fa"), "-o", str(tmp_path / "s")],
                       capture_output=True, text=True)
    assert r.returncode == 1 and "ERROR: input" in r.stderr and "is not a valid BAM" in r.stderr


def test_missing_image_fails(data, tmp_path):
    """Check that a fgbio image that is not loaded fails with an ERROR: message."""
    r = subprocess.run(["bash", str(SCRIPT), "-I", "example.invalid/none:0", "-V", FGBIO_VERSION, "-i", str(data / "in.bam"),
                        "-r", str(data / "ref.fa"), "-o", str(tmp_path / "s")], capture_output=True, text=True)
    assert r.returncode == 1 and "ERROR: docker image" in r.stderr and "is not loaded" in r.stderr


def test_wrong_fgbio_version_fails(data, tmp_path):
    """Check that an image whose fgbio version is not the one asked for fails before any output is written."""
    r = subprocess.run(["bash", str(SCRIPT), "-I", IMAGE, "-V", "4.1.2", "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"),
                        "-o", str(tmp_path / "s"), "-w", str(tmp_path / "work")], capture_output=True, text=True)
    assert r.returncode == 1, r.stderr
    assert f"ERROR: docker image {IMAGE} has fgbio version '{FGBIO_VERSION}', expected 4.1.2" in r.stderr
    assert not list(tmp_path.glob("s.*")) and not (tmp_path / "work").exists()


@pytest.mark.parametrize("version", ["", "latest", "4", "4.1.1-SNAPSHOT"])
def test_malformed_or_missing_fgbio_version_fails(data, tmp_path, version):
    """Check that a missing or malformed -V value fails with an ERROR: message before the image is run."""
    r = subprocess.run(["bash", str(SCRIPT), "-I", IMAGE, "-V", version, "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"),
                        "-o", str(tmp_path / "s")], capture_output=True, text=True)
    assert r.returncode == 1 and "ERROR: " in r.stderr, r.stderr
    assert ("are required" in r.stderr) if not version else ("is not a version number" in r.stderr)
