#!/usr/bin/env python
"""Experiment group for the HCP pipeline."""

import pprint
import os
import re

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.ExperimentGroup import ExpGroup
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import dir_funcs


class HCPGroup(ExpGroup.ExpGroup):
    """Inherits from the ExpGroup object. Some of these functions will be very
    similar to those that are written for the samplesheet object. This is
    unavoidable, unfortunately."""

    def setup(self, args):
        """Do a detailed check of the arguments passed to this function. We
        want to validate the output directory, working directory, BIDS
        directory, and the columns that were passed."""
        valid_args = self._validate(args)
        # This is an ugly hack, but the first column of this list should be the
        # unique key for each *unit*.
        self.columns.extend(['SubjectID,SessionID', 'Group'])
        # Append the extra columns to the default ones
        self.columns.extend(valid_args['extra_column'])
        # The format of these keys will be subject,session: if the session is
        # empty, then it will just be subject,. This is because the full sample
        # unit is a session from a subject.
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
        cannot be found. This does a very CURSORY check that a directory
        looks like a BIDS directory - is there a subject directory?"""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.HCP_BAD_BIDS)
        # Check that there are subjects in the directory, too. These are
        # called sub-[something] and are directories
        sub_pat = re.compile(r'^sub-[a-zA-Z0-9]+$')
        has_subject = False
        for f in contents:
            if re.match(sub_pat, f):
                if dir_funcs.dir_exists(os.path.join(d, f), self.group_logger):
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
        # Regex for the subject directory
        sub_pat = re.compile(r'^sub-[a-zA-Z0-9]+$')
        # regex for the session directory. This is optional if there is just
        # one session.
        ses_pat = re.compile(r'^ses-[a-z-A-Z0-9]+$')
        cont = os.listdir(d)
        sd = {}
        for f in cont:
            if re.match(sub_pat, f):
                # The subject ID is the second part after the '-' character
                subid = f.split('-')[1]
                self.group_logger.debug('Found subject %s', subid)
                # Next, get the session, if it exists
                sub_cont = os.listdir(os.path.join(d, f))
                sessions = []
                for g in sub_cont:
                    if re.match(ses_pat, g):
                        sesid = g.split('-')[1]
                        self.group_logger.debug(
                            'Found session %s for subject %s', sesid, subid)
                        sessions.append(sesid)
                # If sessions is an empty list, then the subject only has one
                # scanning session. Keep the element an empty string.
                if sessions == []:
                    sessions = ['']
                for s in sessions:
                    sd[','.join([subid, s])] = {}
        return sd
