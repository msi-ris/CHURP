#!/usr/bin/env python
"""Define a sub-class of samplesheet that holds the data for the HCP pipeline
"""

import os
import re
import pprint

from GopherPipelines import DieGracefully
from GopherPipelines.SampleSheet import SampleSheet
from GopherPipelines.ArgHandling import set_verbosity


class HCPSampleSheet(SampleSheet.Samplesheet):
    """A sub-class of SampleSheet that holds information for the HCP
    pipelines"""

    def __init__(self, args):
        """Initialize the HCP samplesheet."""
        # Set up a logger
        self.sheet_logger = set_verbosity.verb(args['verbosity'], __name__)
        # Set the column order to be the columns of the sample sheet. This will
        # eventually become the header of the sheet.
        self.column_order.extend([
            'Group',
            'BIDSDir',
            'OutDir',
            'Stages'])
        self.samples = self._parse_bids(args['input_dir'])
        self.stages = ','.join(args['stages'])
        self._resolve_options()
        # Set the sample group memberships based on the expr_groups argument
        self._set_groups(args['expr_groups'])
        return

    def _parse_bids(self, bd):
        """Crawl a BIDS directory structure and return a list of tuples that
        represents the (subject, session) structure inside it."""
        # Regex for the subject directory
        sub_pat = re.compile(r'^sub-[a-zA-Z0-9]+$')
        # regex for the session directory. This is optional if there is just
        # one session.
        ses_pat = re.compile(r'^ses-[a-z-A-Z0-9]+$')
        cont = os.listdir(bd)
        sd = {}
        for f in cont:
            if re.match(sub_pat, f):
                # The subject ID is the second part after the '-' character
                subid = f.split('-')[1]
                self.sheet_logger.debug('Found subject %s', subid)
                # Next, get the session, if it exists
                sub_cont = os.listdir(os.path.join(bd, f))
                sessions = []
                for g in sub_cont:
                    if re.match(ses_pat, g):
                        sesid = g.split('-')[1]
                        self.sheet_logger.debug(
                            'Found session %s for subject %s', sesid, subid)
                        sessions.append(sesid)
                # If sessions is an empty list, then there is only one session
                # for this subject. We make it an empty element
                if sessions == []:
                    sessions = ['']
                for s in sessions:
                    sd[','.join([subid, s])] = {}
        return sd

    def _set_groups(self, groups):
        """If the groups argument is NoneType, then there were no experimental
        groups passed to the pipeline script, and we fill in 'NULL' for each
        group. If it was passed, then we parse it for group memberships. We set
        any overlapping sample groups to be the same value. We report any
        non-overlapping samples as warnings, and complete non-overlap as an
        error."""
        if not groups:
            self.sheet_logger.debug(
                'No groups file passed. All samples are NULL group.')
            for s in self.samples:
                self.samples[s]['Group'] = 'NULL'
            return
        else:
            self.sheet_logger.debug('Parsing %s for groups.', groups)
            csv_gps = {}
            with open(groups, 'r') as f:
                for index, line in enumerate(f):
                    if index == 0:
                        continue
                    elif line.strip() == '':
                        self.sheet_logger.warn(
                            'Line %d in groups csv is empty, skipping.', index)
                    elif len(line.strip().split(',')) < 3:
                        # The subject and session ID are the first two fields,
                        # which uniquely key the "sample"
                        self.sheet_logger.warn(
                            ('Line %d in groups csv has fewer than 3 fields, '
                             'skipping.'), index)
                    else:
                        tmp = line.strip().split(',')
                        csv_gps[tmp[0] + ',' + tmp[1]] = tmp[2]
            self.sheet_logger.debug(
                'CSV experimental groups:\n%s',
                pprint.pformat(csv_gps))
            # Calculate the overlaps
            bids_samples = set(self.samples)
            csv_samples = set(csv_gps)
            bids_only = bids_samples - csv_samples
            csv_only = csv_samples - bids_samples
            overlap = bids_samples & csv_samples
            # Thow an error if there are no overlapping samples
            if len(overlap) == 0:
                self.sheet_logger.error(
                    'The BIDS directory and CSV file appear to mismatch.')
                DieGracefully.die_gracefully(DieGracefully.HCP_NO_SAMP_GPS)
            # Drop warnings if there are exclusive samples
            if len(bids_only) > 0:
                self.sheet_logger.warning(
                    'Samples found that do not have CSV entries, setting to '
                    'NULL group:\n%s',
                    '\n'.join(list(bids_only)))
            if len(csv_only) > 0:
                self.sheet_logger.warning(
                    'Ignoring samples with groups but no data in the BIDS '
                    'directory:\n%s',
                    '\n'.join(list(csv_only)))
            # Iterate through the sample dictionary and set groups
            for s in self.samples:
                gp = csv_gps.get(s, 'NULL')
                self.sheet_logger.debug('Sample %s gets group %s.', s, gp)
                self.samples[s]['Group'] = gp
            return

    def compile(self, bd, od, wd):
        """Iterate through the dictionary of samples, and fill-in the values
        that were supplied by the user."""
        for s in self.samples:
            self.final_sheet[s] = {
                'Group': self.samples[s]['Group'],
                'BIDSDir': str(bd),
                'OutDir': str(od),
                'Stages': self.stages}
        self.sheet_logger.debug(
            'Samplesheet:\n%s',
            pprint.pformat(self.final_sheet))
        return
