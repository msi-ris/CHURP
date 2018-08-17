#!/usr/bin/env python
"""Define a sub-class of samplesheet that holds the data for bulk RNAseq
samples."""

import os
import glob
import re
import pprint
from GopherPipelines.SampleSheet import SampleSheet
from GopherPipelines.ArgHandling import set_verbosity


class BulkRNASeqSampleSheet(SampleSheet.Samplesheet):
    """A sub-class of SampleSheet that holds information for the bulk RNAseq
    pipelines."""

    def __init__(self, args):
        """Initialize the bulk RNAseq samplesheet."""
        # Set up a logger
        self.sheet_logger = set_verbosity.verb(args['verbosity'], __name__)
        # This pipeline takes options for trimmomatic and hisat2
        self.programs.extend(['trimmomatic', 'hisat2'])
        # Set the default trimmomatic options here. This is from YZ's scripts
        self.defaultopts['trimmomatic'] = 'ILLUMINACLIP:' + args['adapters'] + ':4:15:7:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18'
        self.defaultopts['hisat2'] = ''
        # Set the user options here
        self.useropts['trimmomatic'] = args['trimmomatic']
        self.useropts['hisat2'] = args['hisat2']
        if args['no_trim']:
            self.useropts['trim'] = '0'
        else:
            self.useropts['trim'] = '1'
        if args['stranded']:
            self.useropts['stranded'] = '1'
        else:
            self.useropts['stranded'] = '0'
        self.useropts['gtf'] = args['gtf']
        self.useropts['hisat2_idx'] = args['hisat2_idx']
        # Set the column order to be the columns of the sample sheet. This will
        # eventually become the header of the sheet.
        self.column_order.extend([
            'FastqR1files',
            'FastqR2file',
            'OutputDir',
            'WorkingDir',
            'TRIM',
            'trimmomaticOpts',
            'Hisat2index',
            'Hisat2Options',
            'Stranded',
            'AnnotationGTF'])
        self._get_fq_paths(args['fq_folder'])
        self._resolve_options()
        return

    def _get_fq_paths(self, d):
        """Read through the contents of a FASTQ directory and try to build a
        list of samples from it."""
        # Get all files that look like fastq files
        fq_pat = re.compile(r'^.+((.fq(.gz)?$)|(.fastq(.gz)?$))')
        cont = os.listdir(d)
        fqs = []
        for f in cont:
            if re.match(fq_pat, f):
                fqs.append(f)
        self.sheet_logger.debug('Found %i fastq files', len(cont))
        # Then, build the sample names from the fastq names
        snames = set()
        for f in fqs:
            # Strip off the _R1_001.fastq.gz and _R2_001.fastq.gz
            snames.add('_'.join(f.split('_')[:-2]))
        # Cast from set to list
        snames = list(snames)
        self.sheet_logger.debug('Samples found: %s', snames)
        # Now iterate through the list of sample names and populate the sample
        # dictionary
        for s in snames:
            # Use globbing to try to get the full R1 and R2 paths
            r1 = glob.glob(os.path.join(d, s + '_R1_*.f*'))
            r2 = glob.glob(os.path.join(d, s + '_R2_*.f*'))
            # Print warnings if there are multiple files found
            if len(r1) > 1:
                self.sheet_logger.warning(
                    'Sample %s appears to have multiple R1 files associated with it: %s',
                    s, r1)
            if len(r2) > 1:
                self.sheet_logger.warning(
                    'Sample %s appears to have multiple R2 files associated with it: %s',
                    s, r2)
            # If R2, then the sample is paired
            if r2:
                pe = True
            else:
                pe = False
            # Put them into the dictionary
            self.samples[s] = {}
            self.samples[s]['R1'] = r1.pop()
            if pe:
                self.samples[s]['R2'] = r2.pop()
            else:
                self.samples[s]['R2'] = ''
        # Print what we found
        self.sheet_logger.debug('Found samples: %s', pprint.pformat(self.samples))
        return

    def compile(self, od, wd):
        """Iterate through the dictionary of samples, and fill-in the values
        that were supplied by the user."""
        # For each sample...
        for s in self.samples:
            # Key in the values based on the column names
            self.final_sheet[s] = {
                'FastqR1files': self.samples[s]['R1'],
                'FastqR2file': self.samples[s]['R2'],
                'OutputDir': od,
                'WorkingDir': wd,
                'TRIM': self.useropts['trim'],
                'trimmomaticOpts': self.finalopts['trimmomatic'],
                'Hisat2index': self.useropts['hisat2_idx'],
                'Hisat2Options': self.finalopts['hisat2'],
                'Stranded': self.useropts['stranded'],
                'AnnotationGTF': self.useropts['gtf']
                }
        self.sheet_logger.debug('Samplesheet:\n%s', pprint.pformat(self.final_sheet))
        return