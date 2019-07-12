#!/usr/bin/env python
"""Experiment group for the HCP pipeline."""

import pprint
import os

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.ExperimentGroup import ExpGroup
from GopherPipelines.ArgHandling import set_verbosity


class HCPGroup(ExpGroup.ExpGroup):
    """Inherits from the ExpGroup object. Some of these functions will be very
    similar to those that are written for the samplesheet object. This is
    unavoidable, unfortunately."""

    self.columns.extend(['SubjectID', 'SessionID', 'Group'])

    def setup(self, args):
        """Do a detailed check of the arguments passed to this function. We
        want to validate the output directory, working directory, BIDS
        directory, and the columns that were passed."""
        valid_args = self._validate(args)
        # Append the extra columns to the default ones
        self.columns.extend(valid_args['extra_column'])
        self.samples = self._get_bids_data(valid_args['input_dir'])
        self._build_groups()
        return

    def _validate(self, a):
        """Validate the arguments. We want to make sure that the BIDS
        directory is not empty, the columns do not collide with each other, and
        that the names do not have any commas in them."""
        self._validate_bids(a['input_dir'])
        # Drop a warning that specifying extra columns means that there will be
        # some more specialized statistical analysis required
        if a['extra_column']:
            self.group_logger.warning(
                'Specifying additional columns for experimental conditions '
                'is an advanced feature, and will require you to write custom '
                'scripts for statistical analysis. CHURP will do tests on the '
                '"Group" column (present by default), but will not account for'
                ' additional experimental details in your design. This is not '
                'an error message.')
        # Check the experimental columns - first make sure that the names are
        # not duplicated
        tot_col = self.columns + a['extra_column']
        if len(tot_col) != len(set(tot_col)):
            self.group_logger.warning(
                'Duplicate columns specified. This will not cause an error ' +
                'in the Python script, but it may cause an error in any ' +
                'downstream statistical analysis.')
        # Check the supplied columns for bad values
        for e in a['extra_column']:
            if ',' in e:
                self.group_logger.error('Column names cannot contain commas.')
                DieGracefully.die_gracefully(DieGracefully.GROUP_BAD_COL)
        # Turn relative paths into absolute paths
        a['input_dir'] = os.path.realpath(
            os.path.expanduser(a['input_dir']))
        return a

    def _validate_bids(self, d):
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

    def _get_bids_data(self, d):
        """Read the contents of the supplied FASTQ directory and parse out the
        sample names."""
        # regular expression to slice out the samplename from the read name
        samp_re = re.compile(
            r'(_S[0-9]+)?(_L00[1-8])?'
            r'(_R(1|2))?'
            r'_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        # regular expression to get files that look like not-R2 FASTQ files
        fq_re = re.compile(
            r'^.+[^_R2]_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        cont = os.listdir(d)
        sd = {}
        for f in cont:
            if re.match(fq_re, f):
                sn = re.sub(samp_re, '', f)
                sd[sn] = {}
                self.group_logger.debug('Found sample %s', sn)
        return sd
