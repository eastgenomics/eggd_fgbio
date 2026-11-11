# Changelog

## v1.0.0

- First release: fgbio 4.1.1 `ClipBam -c Soft --clip-overlapping-reads true` on a BAM (query-name sort, clip,
  coordinate sort, index), for the AMBER input of the Atlas CNV workflow (DI-4247).
- fgbio image bundled by file ID and checked by SHA-256; samtools from `htslib_suite_asset` 1.24.
- Outputs the clipped BAM and index, fgbio's clipping metrics and a provenance TSV; fails if the record count changes.
- CI: ShellCheck and pytest on a synthetic BAM.
