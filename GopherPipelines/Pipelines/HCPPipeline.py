#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for HCP pipeline"""

import pprint
import os
import glob
import subprocess
import re
import datetime
import getpass

import GopherPipelines
from GopherPipelines import DieGracefully
from GopherPipelines.Pipelines import Pipeline
from GopherPipelines.ArgHandling import set_verbosity
from GopherPipelines.FileOps import default_files
from GopherPipelines.FileOps import dir_funcs
