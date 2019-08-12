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
        self.real_bids = valid_args['input_dir']
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
        a['outdir'] = os.path.realpath(
            os.path.expanduser(str(a['outdir'])))
        a['workdir'] = os.path.realpath(
            os.path.expanduser(str(a['workdir'])))
        # Validoate the stages
        # a['stages'] = self._resolve_stages()
        self.pipe_logger.debug('Input directory: %s', a['input_dir'])
        self._check_input(a['input_dir'])
        return a

    def _resolve_stages(self):
        """Resolve dependencies among the stages. This is complicated - we will
        sort this out later."""
        pass

    def _check_input(self, d):
        """Return an error if the input BIDS directory is not readable or
        cannot be found."""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.HCP_BAD_BIDS)
        # Check that there are subjects in the directory, too. These are
        # called sub-[something] and are directories
        sub_pat = re.compile(r'^sub-.+$')
        has_subject = False
        for f in contents:
            if re.match(sub_pat, f):
                if dir_funcs.dir_exists(os.path.join(d, f), self.pipe_logger):
                    has_subject = True
                    break
        if has_subject:
            return
        else:
            DieGracefully.die_gracefully(DieGracefully.NO_SUBJECTS)
        return

    def _prepare_samplesheet(self):
        """Call the samplesheet build method here. This will build the
        dictionary that will hold all samplesheet data, and then write it into
        the output directory."""
        self._run_checks()
        self.sheet.compile(self.real_bids, self.real_out, self.real_work)
        ss_path = self.sheet.write_sheet(self.real_out, self.pipe_name, '|')
        return ss_path

    def qsub(self):
        """Write the qsub command. We will need the path to the samplesheet,
        the number of samples in the samplesheet, and the scheduler options
        that were passed as arguments. This is defined in the subclass rather
        than in the main Pipeline class because the exact form of the qsub
        command depends on which pipeline we are running."""
        ss = self._prepare_samplesheet()
        # Make the qsub array key
        keyname = default_files.default_array_key(self.pipe_name)
        keyname = os.path.join(self.real_out, keyname)
        if os.path.isfile(keyname):
            self.pipe_logger.warning(
                'Qsub key file %s exists. Overwriting!', keyname)
        try:
            handle = open(keyname, 'w')
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        # The sheet is sorted in this way before it is written to disk, so it
        # should be safe to sort it this way
        handle.write('Qsub.Index\tSampleName\n')
        for index, samplename in enumerate(sorted(self.sheet.final_sheet)):
            handle.write(str(index+1) + '\t' + samplename + '\n')
        handle.flush()
        handle.close()
        # Make the script filename
        pname = default_files.default_pipeline(self.pipe_name)
        pname = os.path.join(self.real_out, pname)
        if os.path.isfile(pname):
            self.pipe_logger.warning(
                'Submission script %s already exists. Overwriting!', pname)
        try:
            handle = open(pname, 'w')
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        # Write the header of the script
        handle.write('#!/bin/bash\n')
        # Write some meta-data lines
        handle.write('# Generated by CHURP version ' +
                     GopherPipelines.__version__ + '\n')
        handle.write('# Generated at ' + GopherPipelines.NOW + '\n')
        handle.write('set -e\n')
        handle.write('set -u\n')
        handle.write('set -o pipefail\n')
        self.pipe_logger.debug(
            'Number of samples: %i', len(self.sheet.final_sheet))
        self.pipe_logger.debug('Samplesheet: %s', ss)
        # Write command to figure out email address of submitting user
        handle.write('user_name="$(id -u -n)"\n')
        handle.write('user_email="${user_name}@umn.edu"\n')
        # This is the command for aligning and cleaning
        qsub_resources = '"mem=' + str(self.mem) + 'mb'
        qsub_resources += ',nodes=1:ppn=' + str(self.ppn)
        qsub_resources += ',walltime=' + str(self.walltime * 3600) + '"'
        # Set the group string here
        if self.group:
            qsub_group = '-A ' + self.group + ' -W "group_list=' + self.group + '"'
        else:
            qsub_group = ''
        qsub_array = '1'
        if len(self.sheet.final_sheet) > 1:
            qsub_array += '-' + str(len(self.sheet.final_sheet))
        # Write a few variables into the header of the script so they are
        # easy to find
        handle.write('KEYFILE=' + '"' + keyname + '"\n')
        handle.write('QSUB_ARRAY=' + '"' + qsub_array + '"\n')
        handle.write('OUTDIR=' + '"' + str(self.real_out) + '"\n')
        handle.write('WORKDIR=' + '"' + str(self.real_work) + '"\n')
        handle.write('SAMPLESHEET=' + '"' + ss + '"\n')
        handle.write('PURGE=' + '"' + self.purge + '"\n')
        handle.write('PIPE_SCRIPT="$(cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/$(basename $0)"\n')
        # Write some echo statements for users' information
        handle.write('echo "Output and logs will be written to ${OUTDIR}"\n')
        handle.write('echo "Emails will be sent to ${user_email}"\n')
        handle.write('echo "Qsub array to samplename key: ${KEYFILE}"\n')
        handle.flush()
        handle.close()
        # Check if we want to automatically submit the script
        if self.nosubmit:
            qsub_dat = None
        else:
            import subprocess
            qsub_cmd = ['bash', pname]
            qsub_proc = subprocess.Popen(
                qsub_cmd,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            qsub_stdout, qsub_stderr = qsub_proc.communicate()
            qsub_dat = (qsub_stdout, qsub_stderr, qsub_proc)
        return (pname, ss, keyname, qsub_dat)
