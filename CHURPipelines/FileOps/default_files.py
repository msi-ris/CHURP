#!/usr/bin/env python
"""Set default filenames for output files."""

import CHURPipelines


def default_pipeline(pipeline):
    """Return a filename for a default pipeline.sh output file."""
    sh_name = '.'.join([
        CHURPipelines.TODAY, CHURPipelines.UNAME, pipeline, 'pipeline.sh'])
    return sh_name


def default_samplesheet(pipeline):
    """Return a filename for a default samplesheet.txt."""
    ss_name = '.'.join([
        CHURPipelines.TODAY, CHURPipelines.UNAME, pipeline,
        'samplesheet.txt'])
    return ss_name


def default_group_csv(pipeline):
    """Return a filename for a default group CSV file."""
    gt_name = '.'.join([
        CHURPipelines.TODAY, CHURPipelines.UNAME, pipeline, 'groups.csv'])
    return gt_name


def default_array_key(pipeline):
    """Return a filename for a qsub array to samplename key."""
    ak_name = '.'.join([
        CHURPipelines.TODAY, CHURPipelines.UNAME, pipeline,
        'qsub_array.txt'])
    return ak_name
