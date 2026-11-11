# Changelog

## v1.0.0

- First release: fgbio `ClipBam -c Soft --clip-overlapping-reads true` on a BAM (query-name sort, clip,
  coordinate sort, index), for the AMBER input of the Atlas CNV workflow (DI-4247).
- fgbio image and version as inputs (`fgbio_docker`, `fgbio_version`): the job fails if the image's `fgbio --version`
  differs, so a new image needs no app change; samtools from `htslib_suite_asset` 1.24.
- Outputs the clipped BAM and index, fgbio's clipping metrics and a provenance TSV; fails if the record count changes
  or a provenance value cannot be read.
- CI on pull requests and pushes to `main`: pytest on a synthetic BAM with samtools/htslib 1.24, ShellCheck 0.11.0 and Ruff
  (including PEP 257 docstring checks).
