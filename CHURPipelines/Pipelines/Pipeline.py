#!/usr/bin/env python
"""Define a general Pipeline object that we will use to build the actual
analysis pipelines. This class will have all methods that are common to the
analysis pipelines that we will build: setting dependencies, checking that
programs exist and checking sample lists."""

import pprint
import os
import sys

import CHURPipelines
from CHURPipelines.ArgHandling import set_verbosity
from CHURPipelines.FileOps import dir_funcs
from CHURPipelines import DieGracefully


class Pipeline(object):
    """Define the Pipeline object."""

    def __init__(self, args):
        """Initialize the pipeline object with user-supplied inputs. The
        general pipeline attributes that get set here are:
            - Output directory
            - Program names for options setting
            - User options dictionary
            - Default optiond dictionary
            - Final options dictionary
            - Path to single sample PBS script"""
        self.logger = set_verbosity.verb(args['verbosity'], __name__)
        self.logger.debug('Passed args:\n%s', pprint.pformat(args))
        self.outdir = args['outdir']
        self.workdir = args['workdir']
        # Set the command log file
        self.cmd_log = args['cmd_log']
        # These are empty, and will get populated by the sub-class.
        self.single_sample_script = ''
        # Set the scheduler resources and PBS options here
        self.group = args['pbs_group']
        self.ppn = args['ppn']
        self.mem = args['mem']
        self.walltime = args['walltime']
        # Set the purge flag
        if args['purge']:
            self.purge = 'true'
        else:
            self.purge = 'false'
        return

    def _run_checks(self):
        """Run the checks on the directories and options that were passed.
        This happens after initialization of sub-classes from the pipeline."""
        self._check_scheduler()
        self._check_dirs()
        self._check_log()
        return

    def _check_log(self):
        """Check that the command log file can be written to."""
        self.logger.info('Checking command log')
        if not self.cmd_log:
            return
        self.logger.debug('Checking command log %s', self.cmd_log)
        cl_fp = os.path.abspath(os.path.expanduser(self.cmd_log))
        cl_dn = os.path.dirname(cl_fp)
        self.logger.debug('Checking %s', cl_dn)
        # First, check that it exists
        if dir_funcs.dir_exists(cl_dn, self.logger):
            # Then, is it writeable?
            if dir_funcs.dir_writeable(cl_dn, self.logger):
                # All good!
                self.logger.debug('Command log dir %s is valid', cl_dn)
            else:
                self.logger.error(
                    'Command log dir %s cannot be written to!', cl_dn)
                DieGracefully.die_gracefully(DieGracefully.BAD_LOGDIR)
        else:
            self.logger.error(
                'Command log dir %s cannot be written to!', cl_dn)
            DieGracefully.die_gracefully(DieGracefully.BAD_LOGDIR)
        return

    def _check_scheduler(self):
        """Check that the scheduler resource requests make sense. ppn should be
        between 1 and 24; mem should be between 2000 and 60000; walltime should
        be between 2h and 96h; and the queues should be one of the valid queues
        on Mesabi or Mangi."""
        try:
            assert self.ppn >= 1 and self.ppn <= 24
        except AssertionError as e:
            self.logger.error(
                'PPN value of %i is invalid! Please specify between 1 and 24.',
                self.ppn)
            DieGracefully.die_gracefully(DieGracefully.BAD_RESOURCES)
        try:
            assert self.mem >= 1 and self.mem <= 60000
        except AssertionError as e:
            self.logger.error(
                'Mem value of %i is invalid! Specify between 1 and 60000.',
                self.mem)
            DieGracefully.die_gracefully(DieGracefully.BAD_RESOURCES)
        try:
            assert self.walltime >= 1 and self.walltime <= 96
        except AssertionError as e:
            self.logger.error(
                'Walltime value of %i is invalid! Specify between 1 and 96.',
                self.walltime)
            DieGracefully.die_gracefully(DieGracefully.BAD_RESOURCES)
        try:
            assert self.msi_queue in CHURPipelines.QUEUES
        except AssertionError as e:
            self.logger.error(
                'Queue %s is not in the allowed list of queues.',
                self.msi_queue)
            DieGracefully.die_gracefully(DieGracefully.BAD_QUEUE)
        return

    def _check_dirs(self):
        """Check that the directories exist and are readable and writeable.
        This will raise an error if we cannot find the fastq directory, or the
        output directory cannot be written to."""
        self.logger.info('Checking directories.')
        # We need to check the output directory and the working directory
        self.logger.debug('Checking output directory %s', self.outdir)
        # First,check that it exists
        if dir_funcs.dir_exists(self.outdir, self.logger):
            # Is it empty?
            if dir_funcs.dir_empty(self.outdir, self.logger):
                # And lastly, is it writeable?
                if dir_funcs.dir_writeable(self.outdir, self.logger):
                    # All good!
                    self.logger.debug('Output dir %s is valid', self.outdir)
                    pass
                else:
                    self.logger.error(
                        'Output dir %s cannot be written to!', self.outdir)
                    DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
            else:
                self.logger.warning(
                    'Output dir %s is not empty! Results may be clobbered.',
                    self.outdir)
        else:
            self.logger.warning(
                'Output dir %s does not exist, making it', self.outdir)
            s = dir_funcs.make_dir(self.outdir, self.logger)
            if not s:
                DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)

        # And do the same for the work directory
        self.logger.debug('Checking working directory %s', self.workdir)
        # First,check that it exists
        if dir_funcs.dir_exists(self.workdir, self.logger):
            # Is it empty?
            if dir_funcs.dir_empty(self.workdir, self.logger):
                # And lastly, is it writeable?
                if dir_funcs.dir_writeable(self.workdir, self.logger):
                    # All good!
                    self.logger.debug('Working dir %s is valid', self.workdir)
                    pass
                else:
                    self.logger.error(
                        'Working dir %s cannot be written to!', self.workdir)
                    DieGracefully.die_gracefully(DieGracefully.BAD_WORKDIR)
            else:
                self.logger.warning(
                    'Working dir %s is not empty!', self.workdir)
        else:
            self.logger.warning(
                'Working dir %s does not exist, making it', self.workdir)
            s = dir_funcs.make_dir(self.workdir, self.logger)
            if not s:
                DieGracefully.die_gracefully(DieGracefully.BAD_WORKDIR)
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
