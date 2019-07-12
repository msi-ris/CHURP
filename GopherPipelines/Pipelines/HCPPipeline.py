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
from GopherPipelines.SampleSheet import HCPSampleSheet
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
        self.sheet = HCPSampleSheet.HCPSampleSheet(valid_args)
        return

    def _validate_args(self, a):
        """Validate arguments for the HCPPipeline object."""
        if not a['input_dir']:
            DieGracefully.die_gracefully(DieGracefully.HCP_INC_ARGS)
        # Convert all of the paths into absolute paths
        a['input_dir'] = os.path.realpath(
            os.path.expanduser(str(a['input_dir'])))
        self.pipe_logger.debug('Input directory: %s', a['input_dir'])
        self._check_input(a['input_dir'])
        return a

    def _check_input(self, d):
        """Return an error if the input BIDS directory is not readable or
        cannot be found."""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.HCP_BAD_BIDS)
        # Check that there is a 'participants.tsv' file
        try:
            open(os.path.join(d, 'participants.tsv'), 'r').close()
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.NO_PARTICIPANTS)
        # Check that there are subjects in the directory, too. These are
        # called sub-[something] and are directories
        sub_pat = re.compile(r'^sub-.+$')
        has_subject = False
        for f in contents:
            if re.match(sub_pat, f, self.pipe_logger):
                if dir_funcs.dir_exists(os.path.join(d, f)):
                    has_subject = True
                    break
        if has_subject:
            return
        else:
            DieGracefully.die_gracefully(DieGracefully.NO_SUBJECTS)
        return
