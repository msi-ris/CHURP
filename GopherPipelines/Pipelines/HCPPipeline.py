#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for HCP pipeline"""

import pprint
import os
import glob
import subprocess
import re
import datetime
import getpass

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import default_files
from GopherPipelines.FileOps import dir_funcs


class HCPPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for HCP pipeline."""
    pipe_name = 'hcp_pipeline'

    def setup(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        # Set up the verbosity and logging
        self.pipe_logger = set_verbosity.verb(args['verbosity'], __name__)
        # First validate the arguments. We do this in the subclass because the
        # dependencies of the arguments are pipeline-specific.
        valid_args = self._validate_args(args)
        self.pipe_logger.debug('New HCPPipeline instance.')
        self.pipe_logger.debug(
            'Validated args:\n%s', pprint.pformat(valid_args))
        self.real_out = valid_args['outdir']
        self.real_work = valid_args['workdir']
        self.nosubmit = valid_args['no_auto_submit']

        # And make a sample sheet from the args
        self.sheet = BulkRNASeqSampleSheet.BulkRNASeqSampleSheet(valid_args)
        return

    def _validate_args(self, a):
        """Validate arguments for the HCPPipeline object."""
        if not a['input_dir']:
            DieGracefully.die_gracefully(DieGracefully.HPC_INC_ARGS)
        # Convert all of the paths into absolute paths
        a['input_dir'] = os.path.realpath(
            os.path.expanduser(str(a['input_dir'])))
        self.pipe_logger.debug('Input directory: %s', a['input_dir'])
        # Check that the adapters and GTF file exist
        try:
            handle = open(a['gtf'], 'r')
            handle.close()
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_GTF)
        if not a['adapters']:
            a['adapters'] = '$TRIMMOMATIC/adapters/all_illumina_adapters.fa'
        else:
            try:
                a['adapters'] = os.path.realpath(
                    os.path.expanduser(str(a['adapters'])))
                handle = open(a['adapters'], 'r')
                handle.close()
            except OSError:
                DieGracefully.die_gracefully(DieGracefully.BAD_ADAPT)
        if a['expr_groups']:
            try:
                handle = open(a['expr_groups'], 'r')
                handle.close()
            except OSError:
                DieGracefully.die_gracefully(DieGracefully.BRNASEQ_BAD_GPS)
        # Validate the FASTQ folder
        self._validate_fastq_folder(a['fq_folder'])
        # Validate the hisat2 index
        self._validate_hisat_idx(a['hisat2_idx'])
        # Sanitize the hisat2 index path
        a['hisat2_idx'] = dir_funcs.sanitize_path(
            a['hisat2_idx'], self.pipe_logger)
        return a
