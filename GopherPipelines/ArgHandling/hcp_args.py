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
        '--output',
        '-o',
        metavar='<Output HCP Dir>',
        dest='output_dir',
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
