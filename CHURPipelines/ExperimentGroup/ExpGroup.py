#!/usr/bin/env python
"""A class to handle the generation of experimental groups for each sample.
This may seem like overkill for now, but eventually we may want to generate
complicated metadata sheets or do some more involved statistical analyses on
groups of samples."""

import sys
import re
import os
import pprint

import CHURPipelines
from CHURPipelines import DieGracefully
from CHURPipelines.ArgHandling import set_verbosity
from CHURPipelines.FileOps import dir_funcs


class ExpGroup(object):
    """The class that will hold the experimetnal group data, as well as write
    the sheet to the output directory."""

    def __init__(self, args):
        """Initialize the object."""
        # Validate the arguments
        self._group_help(args)
        self.cmd_log = args['cmd_log']
        self.group_logger = set_verbosity.verb(args['verbosity'], __name__)
        self.samples = {}
        self.columns = ['SampleName', 'Group']
        self.dest = os.path.realpath(os.path.expanduser(args['outfile']))
        self._prepare_output()
        return

    def _check_log(self):
        """Check that the command log file can be written to."""
        self.group_logger.info('Checking command log')
        if not self.cmd_log:
            return
        self.group_logger.debug('Checking command log %s', self.cmd_log)
        cl_fp = os.path.abspath(os.path.expanduser(self.cmd_log))
        cl_dn = os.path.dirname(cl_fp)
        self.group_logger.debug('Checking %s', cl_dn)
        # First, check that it exists
        if dir_funcs.dir_exists(cl_dn, self.group_logger):
            # Then, is it writeable?
            if dir_funcs.dir_writeable(cl_dn, self.group_logger):
                # All good!
                self.group_logger.debug('Command log dir %s is valid', cl_dn)
            else:
                self.group_logger.error(
                    'Command log dir %s cannot be written to!', cl_dn)
                DieGracefully.die_gracefully(DieGracefully.BAD_LOGDIR)
        else:
            self.group_logger.error(
                'Command log dir %s cannot be written to!', cl_dn)
            DieGracefully.die_gracefully(DieGracefully.BAD_LOGDIR)
        return

    def _group_help(self, args):
        """Do a simple check to see if the pipe_group argument is set. If it
        isn't, print a nice help message."""
        if not args['pipe_group']:
            DieGracefully.die_gracefully(DieGracefully.GROUP_NO_PIPE)
        return

    def _prepare_output(self):
        """Make sure the parent directory of the output directory exists and
        can be written into. Make directories that we have permission to
        make."""
        # Get the dirname of the output file to make sure that the directory
        # exists and can be written to
        par = os.path.dirname(self.dest)
        self.group_logger.info('Checking directories.')
        # We need to check the output directory and the working directory
        self.group_logger.debug('Checking output directory %s', par)
        # First,check that it exists
        if dir_funcs.dir_exists(par, self.group_logger):
            # And lastly, is it writeable?
            if dir_funcs.dir_writeable(par, self.group_logger):
                # All good!
                self.group_logger.debug(
                    'Output dir %s is valid', par)
            else:
                self.group_logger.error(
                    'Output dir %s cannot be written to!', par)
                DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        else:
            self.group_logger.warning(
                'Output dir %s does not exist, making it', par)
            s = dir_funcs.make_dir(par, self.group_logger)
            if not s:
                DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        return

    def _validate_fastq_folder(self, d):
        """Raise an error if the FASTQ directory does not exist or does not
        have any FASTQ files."""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_FASTQ)
        # Check if there is at least one file ending in a standard fastq suffix
        fq_pat = re.compile(r'^.+((.fq(.gz)?$)|(.fastq(.gz)?$))')
        has_fastq = False
        for f in contents:
            if re.match(fq_pat, f):
                has_fastq = True
                break
        if has_fastq:
            return
        DieGracefully.die_gracefully(DieGracefully.EMPTY_FASTQ)
        return

    def _get_sample_names(self, d):
        """Read the contents of the supplied FASTQ directory and parse out the
        sample names."""
        # regular expression to slice out the samplename from the read name
        # This has been updated to also include support for a 4+ nucleotide
        # barcode in the filename
        samp_re = re.compile(
            r'(_S[0-9]+)?'
            r'(_[ATCG]{4,})?'
            r'(_L00[1-8])?'
            r'(_R(1|2))?_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        # regular expression to get files that look like not-R2 FASTQ files
        # in the standard Illumina format
        fq_re = re.compile(
            r'^.+[^_R2]_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        # These regex lines match SRA-style filenames
        sra_re = re.compile(r'^.+_1\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        sra_samp_re = re.compile(r'_(1|2)\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        cont = os.listdir(d)
        sd = {}
        for f in cont:
            if re.match(fq_re, f):
                sn = re.sub(samp_re, '', f)
                sd[sn] = {}
                self.group_logger.debug('Found sample %s', sn)
            elif re.match(sra_re, f):
                sn = re.sub(sra_samp_re, '', f)
                sd[sn] = {}
                self.group_logger.debug('Found SRA sample %s', sn)
        return sd

    def _build_groups(self, dummy='NULL'):
        """Populate the sample dictionary with dummy values."""
        for sn in self.samples:
            for c in self.columns:
                if c == 'SampleName':
                    self.samples[sn][c] = sn
                else:
                    self.samples[sn][c] = dummy
        self.group_logger.debug(
            'Sample groups:\n%s', pprint.pformat(self.samples))
        return

    def write_sheet(self):
        """Write the sheet to the output file."""
        if os.path.isfile(self.dest):
            self.group_logger.warning(
                'Groups file %s exists! Overwriting!', self.dest)
        fh = open(self.dest, 'wt')
        fh.write(','.join(self.columns) + '\n')
        for sn in sorted(self.samples):
            towrite = ','.join(
                [self.samples[sn][c] for c in self.columns]) + '\n'
            fh.write(towrite)
        fh.flush()
        fh.close()
        return

    def write_cmd_log(self):
        """Write the log of the command that was used to invoke the run of
        CHURP. This is separate from the pipeline.sh and samplesheet.txt files
        that get written."""
        if not self.cmd_log:
            return
        cl_fh = open(self.cmd_log, 'at')
        comments = '# CHURP run on {day} {time}\n'
        comments += '# CHURP version: {version}\n'
        comments += '# CWD: {wd}\n'
        cl_fh.write(
            comments.format(
                day=CHURPipelines.TODAY,
                time=CHURPipelines.HUM_TIMESTAMP,
                version=CHURPipelines.__version__,
                wd=os.getcwd())
            )
        cl_fh.write(' '.join(sys.argv) + '\n')
        cl_fh.flush()
        cl_fh.close()
        return
