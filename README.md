<!-- dx-header -->
# eggd_fgbio (DNAnexus Platform App)

Runs [fgbio](https://github.com/fulcrumgenomics/fgbio) tools on DNAnexus. Version 1.0.0 runs one tool, `ClipBam`, to
soft-clip the part of each FR read pair (primary alignments) where the two reads overlap.

## What does this app do?

When the two reads of a pair overlap, a base in the overlap is sequenced twice from the same DNA fragment. Tools that
count reads at a site, such as AMBER, would count it twice. This app soft-clips the overlap so that each base of the
fragment is covered by one read only (see below for the alignments fgbio clips):

```
samtools sort -n  ->  fgbio ClipBam -c Soft --clip-overlapping-reads true  ->  samtools sort  ->  samtools index
```

fgbio clips overlaps between the primary alignments of FR pairs (mates on the same contig, facing each other);
secondary and supplementary alignments and other orientations are not clipped. Soft clipping keeps every base: no read
is removed, and sequences and qualities are unchanged; fgbio may update tags (for example mate information). In rare
alignments (for example a read whose remaining aligned part would start inside a deletion) fgbio can leave a one-base
overlap or unmap a read, changing its flag; on a 103,775,448-record panel BAM the mapped, properly-paired and
mate-mapped counts were the same before and after. The app fails if the output does not hold exactly as many records as
the input.

## What are the typical use cases for this app?

In the Atlas CNV workflow, the tumour BAM passes through this app before AMBER (DI-4247). On BAMs without UMI
consensus, AMBER otherwise counts both reads of an overlapping pair, which widens the B-allele-frequency spread and
leads PURPLE to false low-purity fits (Confluence DV page 4871323786). COBALT, SAGE and CNVkit keep the unclipped BAM.

## What are the inputs?

| Input | Required | Description |
|---|---|---|
| `bam` | yes | Paired-end BAM. Its index is not needed. |
| `reference_fasta` | yes | FASTA the BAM is aligned to, plain or bgzipped (Atlas: GRCh38 GIABv3, `project-Fkb6Gkj433GVVvj73J7x8KbV:file-Gb73P8Q4XGyX9QZ0g9qp61xV`). |
| `reference_fai` | yes | `<reference_fasta>.fai`. |
| `reference_gzi` | if bgzipped | `<reference_fasta>.gzi`. |
| `output_prefix` | no | Output file prefix; defaults to the BAM name without `.bam`. |

A sequence dictionary for the reference is made on the worker with `samtools dict`.

## What are the outputs?

| Output | File |
|---|---|
| `clipped_bam` | `<prefix>.clipped.bam`, coordinate-sorted |
| `clipped_bai` | `<prefix>.clipped.bam.bai` |
| `clipbam_metrics` | `<prefix>.clipbam_metrics.txt`, fgbio's clipping metrics |
| `provenance` | `<prefix>.provenance.tsv` (columns `field`, `value`, `qc_failed`): app, fgbio image ID and version, samtools version, settings, runtime, record counts, and `samtools flagstat` before (`flagstat_in:`) and after (`flagstat_out:`) with QC-passed and QC-failed counts |

The job fails, before any output is published, if any provenance value cannot be read.

## How to run this app from the command line?

```bash
R=project-Fkb6Gkj433GVVvj73J7x8KbV
dx run eggd_fgbio \
  -ibam=project-xxxx:file-xxxx \
  -ireference_fasta=$R:file-Gb73P8Q4XGyX9QZ0g9qp61xV \
  -ireference_fai=$R:file-Gb735V84XGyYkXYPqGqgz14j \
  -ireference_gzi=$R:file-Gb735Vj4XGyj2ZPbzjfYxZxX
```

The default instance is `mem2_ssd1_v2_x16`; a 6 GB panel BAM (about 100 million reads) takes about 30 minutes.

## Dependencies

- **fgbio 4.1.1**: the image `quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace`
  (tag `4.1.1--hdfd78af_0`), saved with `docker save` and bundled with the app (`runSpec.bundledDepends`). The app checks
  the saved image's SHA-256 (`6740248d…99c6`) before loading it, so no image is pulled at run time.
- **samtools**: `htslib_suite_asset` 1.24 (`record-J9QzZyj487v9Gj2QgP7vK47g`).

## Tests

`tests/test_clip_overlap.py` runs `resources/home/dnanexus/clip_overlap.sh` on a small synthetic BAM (generated in
the test). It checks that overlapping pairs are clipped and no longer overlap; that flags, sequences and qualities are
unchanged for every read, and positions and CIGARs for reads that do not overlap their mate; that the record count is
kept; that a bgzipped reference gives the same records as a plain one; that the provenance has its fields; and that bad
inputs and an unreadable provenance value fail with an `ERROR:` message. It needs samtools, bgzip and the fgbio image in
the local docker daemon:

```bash
docker pull quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace
docker tag quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace \
  quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0
python -m pip install -r requirements-dev.txt && python -m pytest tests
```

CI (`.github/workflows/ci.yml`) runs these tests with samtools/htslib 1.24 (the app's version, built from the release
tarball), ShellCheck 0.11.0 and Ruff on every push and pull request.

## Making the bundled image

```bash
docker pull quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace
docker tag  quay.io/biocontainers/fgbio@sha256:a1b0630a647886d67bab11992370502e6c9239d1c4fe01f4d0125908980e1ace \
  quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0
mkdir -p bundle/home/dnanexus/fgbio_image
docker save quay.io/biocontainers/fgbio:4.1.1--hdfd78af_0 -o bundle/home/dnanexus/fgbio_image/fgbio_4.1.1--hdfd78af_0.docker.tar
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime='2026-11-11 00:00Z' -C bundle -cf - home \
  | gzip -n -6 > fgbio_4.1.1--hdfd78af_0.bundle.tar.gz
```

The bundle is uploaded to DNAnexus and its file ID set in `dxapp.json`; the SHA-256 of the inner `.docker.tar` is set
in `src/code.sh`. The v1.0.0 bundle is `file-JFp2PYj453Kg6KGFB694q6bf`, in 001_Reference at
`/app_assets/fgbio/fgbio_4.1.1/` (a copy into another project keeps the file ID). `dx build --app` copies it into the
app, so the built app does not depend on the source file.

## Release

Following the Bioinformatics Development Manual v3: build the app from the reviewed commit, test it, record the
evidence on the validation report in Confluence (DV), and after sign-off publish that same app
(for v1.0.0: `dx publish eggd_fgbio/1.0.0`), merge the release branch into `main` and tag the release.
