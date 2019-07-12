#!/usr/bin/env python
"""Add arguments for the UMII HCP pipeline."""

import argparse
from GopherPipelines.FileOps import default_dirs


def add_args(ap):
    """Takes an ArgumentParser object, and adds some arguments to it. These
    will be for the bulk RNAseq pipeline. This function returns None; it is
    called for the side-effect of adding arguments to the parser object."""
    ap_req = ap.add_argument_group(
        title='Required Arguments')
    ap_req.add_argument(
        '--input',
        '-i',
        metavar='<Input BIDS Dir>',
        dest='input_dir',
        help='Input data directory, in BIDS structure.')
    # Make an argument group for optional arguments
    ap_opt = ap.add_argument_group(
        'Optional Arguments')
    ap_opt.add_argument(
        '--help',
        '-h',
        help='Show this help message and exit.',
        action='help')
    ap_opt.add_argument(
        '--verbosity',
        '-v',
        metavar='<loglevel>',
        dest='verbosity',
        help=('How much logging output to show. '
              'Choose one of "debug," "info," or "warn." Default: warn'),
        choices=['debug', 'info', 'warn'],
        default='warn')
    ap_opt.add_argument(
        '--working-dir',
        '-d',
        metavar='<Working Dir>',
        dest='workdir',
        help=('Working directory for intermediate files. '
              'Defaults to global scratch')),
    ap_opt.add_argument(
        '--output',
        '-o',
        metavar='<Output HCP Dir>',
        dest='outdir',
        help='Output directory. Defaults to global scratch',
        default=default_dirs.default_output('hcp'))
    ap_opt.add_argument(
        '--stages',
        '-s',
        metavar='<Stages>',
        nargs='*',
        default=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer', 'fMRIVolume',
                 'fMRISurface', 'multiICAFIX', 'MSMAll',
                 'DiffusionPreProcessing'],
        help='Stages to run in the HCP pipeline.',
        choices=['PreFreeSurfer', 'FreeSurfer', 'PostFreeSurfer', 'fMRIVolume',
                 'fMRISurface', 'multiICAFIX', 'MSMAll',
                 'DiffusionPreProcessing'])
    ap_opt.add_argument(
        '--no-submit',
        help=('Do not automatically submit pipeline jobs. Use this when you '
              'want to make changes to the pipeline script or samplesheet '
              'before running.'),
        dest='no_auto_submit',
        action='store_true',
        default=False)
    ap_opt.add_argument(
        '--purge',
        dest='purge',
        help='If supplied, overwrite files from pervious runs.',
        action='store_true',
        default=False)
    # Make an argument for scheduler options
    ap_sched = ap.add_argument_group(
        'Scheduler Options')
    ap_sched.add_argument(
        '--group',
        '-A',
        metavar='<MSI group>',
        dest='pbs_group',
        help='MSI group to charge for SU usage. Defaults to primary group.',
        type=str,
        default=None)
    ap_sched.add_argument(
        '--ppn',
        '-p',
        metavar='<procs per node>',
        dest='ppn',
        help='Processors to allocate for each task. Defaults to 6.',
        type=int,
        default=6)
    ap_sched.add_argument(
        '--mem',
        '-m',
        metavar='<mem per job (Mb)>',
        dest='mem',
        help=('Memory, in megabytes, to allocate to each task. Must be at '
              'least 12000 (12GB). Default: 12000'),
        type=int,
        default=12000)
    ap_sched.add_argument(
        '--walltime',
        '-w',
        metavar='<wall clock time (hours)>',
        dest='walltime',
        help=('Walltime, in hours, to allocate to each task. Must be at least '
              '2 hours. Defaults to 12 hours.'),
        type=int,
        default=12)
    return
