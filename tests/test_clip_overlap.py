"""Tests for resources/home/dnanexus/clip_overlap.sh on a small synthetic BAM.

Needs samtools on PATH and the fgbio image (quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0) in the local docker daemon;
CI pulls it by digest. The reference and reads are generated here (random sequence, fixed seed), so no patient data.
"""
import random
import re
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "resources/home/dnanexus/clip_overlap.sh"
IMAGE = "quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0"
READ_LEN = 100
# name: (pos1, pos2), 1-based leftmost positions of read 1 (forward) and read 2 (reverse) on chr1
PAIRS = {
    "overlap_40": (1001, 1061),       # reads overlap by 40 bases
    "overlap_full": (2001, 2001),     # reads cover the same 100 bases
    "abutting": (3001, 3101),         # no overlap: read 2 starts after read 1 ends
    "apart": (4001, 4301),            # far apart
}


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    rnd = random.Random(4247)
    seq = {c: "".join(rnd.choice("ACGT") for _ in range(6000)) for c in ("chr1", "chr2")}
    fa = d / "ref.fa"
    fa.write_text("".join(f">{c}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n" for c, s in seq.items()))
    run(["samtools", "faidx", str(fa)])
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
    out = {}
    for line in run(["samtools", "view", str(bam)]).stdout.splitlines():
        f = line.split("\t")
        out[(f[0], int(f[1]) & 0xC0)] = (int(f[3]), f[5])
    return out


def ref_end(pos, cigar):
    return pos + sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X") - 1


@pytest.fixture(scope="module")
def clipped(data, tmp_path_factory):
    o = tmp_path_factory.mktemp("out")
    run(["bash", str(SCRIPT), "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"), "-o", str(o / "s"), "-t", "2",
         "-m", "100M", "-x", "1", "-w", str(o / "work")])
    return o


def test_outputs_exist(clipped):
    for suffix in (".clipped.bam", ".clipped.bam.bai", ".clipbam_metrics.txt", ".provenance.tsv"):
        assert (clipped / f"s{suffix}").stat().st_size > 0, suffix


def test_record_count_unchanged(data, clipped):
    assert len(records(clipped / "s.clipped.bam")) == len(records(data / "in.bam")) == 2 * len(PAIRS) + 2


def test_overlapping_pairs_no_longer_overlap(clipped):
    rec = records(clipped / "s.clipped.bam")
    for name in ("overlap_40", "overlap_full"):
        (p1, c1), (p2, c2) = rec[(name, 0x40)], rec[(name, 0x80)]
        assert "S" in c1 + c2, (name, c1, c2)
        assert ref_end(p1, c1) < p2 or ref_end(p2, c2) < p1, (name, p1, c1, p2, c2)
    # 40 overlapping bases are split between the two reads
    (p1, c1), (p2, c2) = rec[("overlap_40", 0x40)], rec[("overlap_40", 0x80)]
    assert ref_end(p1, c1) - p1 + 1 + ref_end(p2, c2) - p2 + 1 == 2 * READ_LEN - 40


def test_other_reads_unchanged(data, clipped):
    before, after = records(data / "in.bam"), records(clipped / "s.clipped.bam")
    for key in before:
        if key[0] not in ("overlap_40", "overlap_full"):
            assert after[key] == before[key], key


def test_provenance(clipped):
    prov = dict(line.split("\t", 1) for line in (clipped / "s.provenance.tsv").read_text().splitlines()[1:])
    assert prov["fgbio_version"].strip() == "4.1.1"
    assert prov["fgbio_image"].strip() == IMAGE
    assert prov["records_in"] == prov["records_out"]


def test_missing_fai_fails(data, tmp_path):
    (tmp_path / "ref.fa").write_text((data / "ref.fa").read_text())
    r = subprocess.run(["bash", str(SCRIPT), "-i", str(data / "in.bam"), "-r", str(tmp_path / "ref.fa"),
                        "-o", str(tmp_path / "s")], capture_output=True, text=True)
    assert r.returncode == 1 and "ERROR: reference index" in r.stderr


def test_not_a_bam_fails(data, tmp_path):
    bad = tmp_path / "bad.bam"
    bad.write_text("not a bam\n")
    r = subprocess.run(["bash", str(SCRIPT), "-i", str(bad), "-r", str(data / "ref.fa"), "-o", str(tmp_path / "s")],
                       capture_output=True, text=True)
    assert r.returncode == 1 and "is not a valid BAM" in r.stderr


def test_missing_image_fails(data, tmp_path):
    r = subprocess.run(["bash", str(SCRIPT), "-i", str(data / "in.bam"), "-r", str(data / "ref.fa"), "-o",
                        str(tmp_path / "s"), "-I", "example.invalid/none:0"], capture_output=True, text=True)
    assert r.returncode == 1 and "is not loaded" in r.stderr
