#!/usr/bin/env python
"""Experiment group for the bulk RNAseq pipeline."""

import pprint
import os
import pandas as pd

import CHURPipelines
from CHURPipelines import DieGracefully
from CHURPipelines.ExperimentGroup import ExpGroup
from CHURPipelines.ArgHandling import set_verbosity


class BulkRNAseqGroup(ExpGroup.ExpGroup):
    """Inherits from the ExpGroup object. Some of these functions will be very
    similar to those that are written for the bulk RNAseq samplesheet object.
    This is unavoidable, unfortunately."""

    def __init__(self, args):
        """Do a detailed check of the arguments passed to this function. We
        want to validate the output directory, working directory, FASTQ
        directory, and the columns that were passed."""
        ExpGroup.ExpGroup.__init__(self, args)
        valid_args = self._validate(args)
        self._ensure_dest_suffix_xlsx()
        self.samples = self._get_sample_names(valid_args['fq_folder'])
        self._build_groups()
        return
    
    def _ensure_dest_suffix_xlsx(self):
        """ If a destination suffix other than xlsx was given, change it """
        dest_parts = self.dest.split('.')
        no_suffix_parts = '.'.join(dest_parts[:-1])
        suffix = dest_parts[-1]
        self.dest = no_suffix_parts + '.xlsx'
        if not suffix in ['xls', 'xlsx']:
            self.group_logger.warning(
                'Groups file must end in xlsx, converting to : %s', self.dest)            
    
    def write_sheet(self):
        """Write a stub excel spreadsheet to the output file."""
        if os.path.isfile(self.dest):
            self.group_logger.warning(
                'Groups file %s exists! Overwriting!', self.dest)
        
        # populate the group sheet with the sample names
        groups = pd.DataFrame({'SampleName': list(self.samples.keys()),
                    'Group': 'NULL'})
        contrasts = pd.DataFrame({'Comparison_Name' : [],
                      'Reference_Group' : [],
                      'Test_Group' : []})
        # save the group and contrast sheets
        with pd.ExcelWriter(self.dest) as writer:  
            groups.to_excel(writer, sheet_name='groups', index = False)
            contrasts.to_excel(writer, sheet_name='contrasts', index = False)
        return

    def _validate(self, a):
        """Validate the arguments. We want to make sure that the FASTQ
        directory is not empty, the columns do not collide with each other, and
        that the names do not have any commas in them."""
        self._validate_fastq_folder(a['fq_folder'])
        # Drop a warning that specifying extra columns means that there will be
        # some more specialized statistical analysis required
        # Check the experimental columns - first make sure that the names are
        # not duplicated
        tot_col = self.columns
        if len(tot_col) != len(set(tot_col)):
            self.group_logger.warning(
                'Duplicate columns specified. This will not cause an error ' +
                'in the Python script, but it may cause an error in any ' +
                'downstream statistical analysis.')
        # Turn relative paths into absolute paths
        a['fq_folder'] = os.path.realpath(
            os.path.expanduser(a['fq_folder']))
        return a
