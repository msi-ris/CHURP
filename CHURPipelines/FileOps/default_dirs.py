#!/usr/bin/env python
"""Set the default directories based on a base path, the user name, and the
date."""

import pathlib
import CHURPipelines

# Set the scratch directory as a constant
SCRATCH = pathlib.Path('/panfs/roc/scratch')


def default_workdir(pipeline):
    """Return a pathlib.Path object that points to a default working directory.
    This will be in global scratch at
    $SCRATCH/username/YYYY-MM-DD.pipeline_id.work"""
    # Make the basename of the directory from YMD, time, and pipeline name
    wk_base = CHURPipelines.TODAY + CHURPipelines.TIMESTAMP + '.' + pipeline + '.work'
    # '/' operator joins pathlib.Path objects
    d = SCRATCH / CHURPipelines.UNAME / wk_base
    return d


def default_output(pipeline):
    """Return a pathlib.Path to a default output scratch directory. This will be
    $SCRATCH/username/YYYY-MM-DD.pipeline_id."""
    outbase = CHURPipelines.TODAY + CHURPipelines.TIMESTAMP + '.' + pipeline
    # Path objects are joined with the / operator
    d = SCRATCH / CHURPipelines.UNAME / outbase
    return d
