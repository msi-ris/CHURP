#!/usr/bin/env python3
"""CHURP control script, for running standard analyses of various types of
high-throughput sequencing data analyses. This script is written to work on the
University of Minnesota Supercomputing Institute clusters with data from the
University of Minnesota Genomics Centre. The following analysis pipelines are
supported:
    - group_template
    - bulk_rnaseq
Questions should be directed to help@msi.umn.edu.
Version: 0.2.2
2021-02-03
"""

# Check the Python version
try:
    import sys
    assert sys.version_info[0] == 3
except AssertionError:
    sys.stderr.write('Error - this pipeline requires Python 3!\n')
    sys.exit(2)


# Import standard library modules here
import sys
import os
import argparse
# Import the package here
try:
    import CHURPipelines
    from CHURPipelines.ArgHandling import args
    from CHURPipelines import DieGracefully
except ImportError:
    sys.stderr.write('Error - CHURPipelines package cannot be not found!\n')
    sys.exit(1)


def sp_dbs(args):
    """This function will read and display a list of available species
    databases that can be used as targets for various pipelines."""
    from CHURPipelines import FavoriteSpecies
    # Print a nice message to describe the table we are showing
    msg = """Favorite Species

The species listed in the table below are common genomics models for which we
have provided convenient shortcuts in CHURP. Use the value in the "Alias" column
to automatically set the relevant genome indices and annotation files in other
CHURP pipelines. These names are *case sensitive.*\n\n"""
    sys.stderr.write(msg)
    # Print a header for the table
    hdr = [
        'Alias '.rjust(15),
        ' Binomial Name'.ljust(28),
        ' Assembly/Annotation'.ljust(33)
        ]
    sys.stderr.write('|' + '|'.join(hdr) + '|\n')
    sys.stderr.write('|' + '-'*15 + '|' + '-'*28 + '|' + '-'*33 + '|\n')
    # Iterate through the FAVE_ASM dict to print a nice table that
    # gives the information for the "favorite species"
    for f_sp, sp_dat in sorted(FavoriteSpecies.FAVE_ASM.items()):
        common_name = f_sp
        sci_name = sp_dat[1]
        anno_ver = sp_dat[2]
        to_print = [
            common_name.rjust(14) + ' ',
            ' ' + sci_name.replace('_', ' ').ljust(27),
            ' ' + anno_ver.ljust(32)
            ]
        sys.stderr.write('|' + '|'.join(to_print) + '|\n')
    return


def expr_group(args):
    """This function will generate experimental group CSV files for input into
    the various pipelines."""

    def brnaseq_group(a):
        """Sub-function for calling the bulk RNAseq group template."""
        from CHURPipelines.ExperimentGroup import BulkRNAseqGroup
        eg = BulkRNAseqGroup.BulkRNAseqGroup(args)
        eg.write_sheet()
        eg.write_cmd_log()
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_GROUP_SUCCESS,
            eg.dest)
        return

    grp_cmd = {
        'bulk_rnaseq': brnaseq_group
        }
    try:
        grp_cmd[args['pipe_group']](args)
    except KeyError:
        # If we get here, then we should get trapped by the no-args help func
        from CHURPipelines.ExperimentGroup import ExpGroup
        e = ExpGroup.ExpGroup(args)
    return


def brnaseq(args):
    """This function loads the bulk RNAseq pipeline module, and runs through
    the steps for bulk RNAseq analysis."""
    from CHURPipelines.Pipelines import BulkRNAseq
    p = BulkRNAseq.BulkRNAseqPipeline(args)
    pipeline_fname, samplesheet_fname, key_name, qsub_dat = p.qsub()
    p.write_cmd_log()
    if not qsub_dat:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUCCESS,
            pipeline_fname,
            samplesheet_fname,
            key_name)
    elif qsub_dat[2].returncode != 0:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUBMIT_FAIL,
            qsub_dat)
    else:
        DieGracefully.die_gracefully(
            DieGracefully.BRNASEQ_SUBMIT_OK,
            pipeline_fname,
            samplesheet_fname,
            key_name,
            qsub_dat)
    return


def main():
    """The main function. This function is a very high-level function, and
    it should really only have the logic and structure of the pipeline that is
    invoked. The modules and accessory functions within the package should
    do all the actual work."""
    if not sys.argv[1:]:
        args.usage()
        sys.exit(3)
    else:
        pipe_args = args.pipeline_args()
        cmd = {
            'bulk_rnaseq': brnaseq,
            'group_template': expr_group,
            'show_faves': sp_dbs
            }
        cmd[pipe_args['pipeline']](pipe_args)
    return


main()
