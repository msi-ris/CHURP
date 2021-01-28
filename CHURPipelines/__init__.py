#!/usr/bin/env python
"""CHURP package for building samplesheets and qsub commands to
run analysis on the Minnesota Supercomputing Institute (MSI) cluster."""

# set dunder (double-underscore) variables
__version__ = '0.2.2-slurm'
__date__ = '2020-10-19'

import datetime
import getpass

# Define today's date and the current time
curr_time = datetime.datetime.now()
TODAY = datetime.date.isoformat(datetime.datetime.today())
# A time stamp to append to the date stamp, in ISO format
TIMESTAMP = 'T{t.hour:02}{t.minute:02}{t.second:02}'.format(t=curr_time)
NOW = curr_time.isoformat(sep=' ', timespec='seconds')
UNAME = getpass.getuser()

# Define the allowable queues
QUEUES = ['small', 'amdsmall']
