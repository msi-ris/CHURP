# CHURP Changelog
2022-05-06

## [TODO]

## [0.2.3] 2022-05-06
Patch-level update for PURR, the CHURP bulk RNAseq workflow. This update adds support for Agate, the new cluster at MSI.

**Note that if you are using Agate, you should request more memory than you normally would for Mangi or Mesabi.** We recommend 4000MB per CPU.

### Added
- `agsmall` (Agate!) and `amd512` partitions to the allowable target list
- `--command-log` option to the `group_template` subcommand
- `--tmp` option to `bulk_rnaseq` to allow users to set a custom local scratch allocation
- Reference genome path and genome annotation path added to the HTML report. This is to further document the reference datasets uses in analysis.
- This changelog!

### Modified
- Job submissions now use `--mail-type=BEGIN,END,FAIL` rather than `ALL`
- Job scripts reset `$PATH` before running any analyses. This will hopefully resolve errors when users have custom `$PATH` definitions that are inherited by the job environment.
- Global scratch path has been updated
- Internal class structure changed slightly; this has no bearing on the user experience

### Removed
- Nothing removed.

## [0.2.2] 2020-04-08
Patch-level update for PURR, the CHURP bulk RNAseq workflow.

### Added
- `-q` flag for specifying which queue to submit jobs. CHURP can now send tasks to the [Mangi](https://www.msi.umn.edu/queues) queues.
- "Large dataset" support; the HTML report now breaks large data sets into
  chunks for display.
- Custom exit codes for each step for easier identification of where PURR fails. See the[documentation](https://github.umn.edu/MSI-RIS/CHURP/wiki/PURR-Manual-Page#exit-codes) for a description of each code.

### Modified
- Strand-specificity and gDNA metrics for single-read data have been fixed.
- HTML report is still produced, even if DEG testing does not take place.
- PURR no longer aborts when RNASeQC cannot be run. This means that non-Ensembl GTFs can be used.
- PURR quits when it detects a mix of single-read and paired-end data in the same FASTQ folder.

### Removed
- Nothing major removed from PURR
