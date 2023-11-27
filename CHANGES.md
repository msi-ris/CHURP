# CHURP Changelog
2023-11-27

## [TODO]

## [0.2.4] 2023-11-27
Patch-level release of CHURP and PURR (RNAseq pipeline). This update adds a new
subcommand and many small tweaks to the output.

# Added
- `genome_aliases` subcommand to CHURP to list aliases for common genomics
  models, see below
- Routines to PURR that color-code plots by experimental group
  (e.g., treatment) and a batching variable
- Downloadable counts matrices to PURR HTML report
- Downloadable pipeline.sh script file to PURR HTML report
- Downloadable samplesheet.txt file to PURR HTML report

# Modified
- PURR HTML report is much smaller and less redundant now
- GTF md5sum is now reported in PURR HTML report
- Uninformative rRNA plots from RNASeQC were removed from PURR HTML report
- `featureCounts` raw counts matrix no longer has transcript-level metadata;
  this makes the resulting file much smaller and easier to store/share
- Automated DEG test outputs from PURR now report all genes, not just those
  with FDR<0.05
- Old package names have been standardized

# Removed
- Nothing else removed

# Bugs Fixed
- #114: DEGs should now have the correct polarity as described in the PURR
  report

## Genome Aliases
When you run `python churp.py genome_aliases`, a table with commonly-used
genomics models is printed. These are genomes that MSI maintains in a local
library ("bioref;" https://www.msi.umn.edu/content/bioref). You can use the
value in the `Alias` column to act as a shorthand for HISAT2 index and GTF
annotation file when using PURR. For example:

```
python churp.py bulk_rnaseq -r corn ...
```

is equivalent to

```
python churp.py bulk_rnaseq \
    -x /common/bioref/ensembl/plants/Zea_mays-56/Zm-B73-REFERENCE-NAM-5.0/hisat2/genome \
    -g /common/bioref/ensembl/plants/Zea_mays-56/Zm-B73-REFERENCE-NAM-5.0/annotation/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf \
    ...
```

The HTML report will still record the full path to the GTF and HISAT2 index,
not just `corn`. Use the `--command-log` option to save the command exactly as
written to a separate file for your records.

## [0.2.3] 2022-05-06
Patch-level update for PURR, the CHURP bulk RNAseq workflow. This update adds
support for Agate, the new cluster at MSI.

**Note that if you are using Agate, you should request more memory than you
normally would for Mangi or Mesabi.** We recommend 4000MB per CPU.

### Added
- `agsmall` (Agate!) and `amd512` partitions to the allowable target list
- `--command-log` option to the `group_template` subcommand
- `--tmp` option to `bulk_rnaseq` to allow users to set a custom local scratch
  allocation
- Reference genome path and genome annotation path added to the HTML report.
  This is to further document the reference datasets uses in analysis.
- This changelog!

### Modified
- Job submissions now use `--mail-type=BEGIN,END,FAIL` rather than `ALL`
- Job scripts reset `$PATH` before running any analyses. This will hopefully
  resolve errors when users have custom `$PATH` definitions that are inherited
  by the job environment.
- Global scratch path has been updated
- Internal class structure changed slightly; this has no bearing on the user
  experience

### Removed
- Nothing removed.

## [0.2.2] 2020-04-08
Patch-level update for PURR, the CHURP bulk RNAseq workflow.

### Added
- `-q` flag for specifying which queue to submit jobs. CHURP can now send tasks
  to the [Mangi](https://www.msi.umn.edu/queues) queues.
- "Large dataset" support; the HTML report now breaks large data sets into
  chunks for display.
- Custom exit codes for each step for easier identification of where PURR
  fails. See the [documentation](https://github.umn.edu/MSI-RIS/CHURP/wiki/PURR-Manual-Page#exit-codes)
  for a description of each code.

### Modified
- Strand-specificity and gDNA metrics for single-read data have been fixed.
- HTML report is still produced, even if DEG testing does not take place.
- PURR no longer aborts when RNASeQC cannot be run. This means that non-Ensembl
  GTFs can be used.
- PURR quits when it detects a mix of single-read and paired-end data in the
  same FASTQ folder.

### Removed
- Nothing major removed from PURR
