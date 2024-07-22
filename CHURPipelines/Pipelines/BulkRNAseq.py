#!/usr/bin/env python
"""Define a sub-class of the Pipeline class for bulk RNAseq analysis."""

import pprint
import os
import glob
import subprocess
import re
import pandas as pd

import CHURPipelines
from CHURPipelines import DieGracefully
from CHURPipelines import FavoriteSpecies
from CHURPipelines.Pipelines import Pipeline
from CHURPipelines.SampleSheet import BulkRNASeqSampleSheet
from CHURPipelines.ArgHandling import set_verbosity
from CHURPipelines.FileOps import default_files
from CHURPipelines.FileOps import dir_funcs


class BulkRNAseqPipeline(Pipeline.Pipeline):
    """Sub-class of Pipeline for bulk RNAseq analysis."""

    # Define the pipeline name here
    pipe_name = 'bulk_rnaseq'

    def __init__(self, args):
        """Initialize the pipeline object. We will call the general
        Pipeline.__init__() here, as well as set some specific pipeline
        attributes."""
        Pipeline.Pipeline.__init__(self, args)
        # Set up the verbosity and logging
        self.pipe_logger = set_verbosity.verb(args['verbosity'], __name__)
        # First validate the arguments. We do this in the subclass because the
        # dependencies of the arguments are pipeline-specific.
        valid_args = self._validate_args(args)
        self.valid_args = valid_args
        self.pipe_logger.debug('New BulkRNAseqPipeline instance.')
        self.pipe_logger.debug(
            'Validated args:\n%s', pprint.pformat(valid_args))
        self.real_out = valid_args['outdir']
        self.real_work = valid_args['workdir']
        self.nosubmit = valid_args['no_auto_submit']
        # Set the boolean summary_only option
        if valid_args['summary_only']:
            self.summary_only = 'true'
        else:
            self.summary_only = 'false'

        # Set the minimum gene length
        self.min_gene_len = str(valid_args['mingene'])
        # And the minimum depth
        self.min_cts = str(valid_args['mincts'])
        # Set the subsampling level
        self.rrna_screen = str(valid_args['rrna_screen'])
        self.subsample = str(valid_args['subsample'])
        # Set the destination queue
        self.msi_queue = str(valid_args['msi_queue'])

        # And make a sample sheet from the args
        self.sheet = BulkRNASeqSampleSheet.BulkRNASeqSampleSheet(valid_args)

        # Set the paths to the single sample analysis script. This will be
        # submitted to the scheduler. This is a little ugly, but because the
        # package contains the shell script data, it should be OK to define a
        # relative path. __file__ is the currently running script,
        # os.path.realpath() gives the full path, with symlinks resolved. We
        # then rsplit() the string returned by realpath() to get the base dir
        # of the CHURP scipt. Woof.
        self.single_sample_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'bulk_rnaseq_single_sample.sh')
        self.summary_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'PBS',
            'run_summary_stats.sh')
        self.de_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'R_Scripts',
            'summarize_bulk_rnaseq.R')
        self.report_script = os.path.join(
            os.path.realpath(__file__).rsplit(os.path.sep, 3)[0],
            'R_Scripts',
            'bulk_rnaseq_report.Rmd')
        return

    def _validate_args(self, a):
        """Validate arguments for the BulkRNAseqPipeline object. We define it
        in this file because it only really needs to be accessible to this
        subclass. Argument dependencies are pipeline-specific. For the bulk
        RNAseq analysis pipeline:
            - FASTQ dir and UMGC sheet are mutually exclusive
            - Organism and HISAT2 index + GTF are mutually exclusive
            - Check that the HISAT2 index is complete
        Further, sanitize the paths of the output dir, working dir, and hisat2
        index.
        """
        # Check the completeness of the argument dictionary. -f must be
        # specified and (-x and -g) or -r must be specified. After checking the
        # FASTQ folder, check the helper commands
        if not a['fq_folder']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif not ((a['hisat2_idx'] and a['gtf']) or a['organism']):
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_INC_ARGS)
        elif (a['hisat2_idx'] or a['gtf']) and a['organism']:
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_CONFLICT)
        # Set the hisat index and gtf if the 'organism' option was supplied
        if a['organism']:
            try:
                org = a['organism']
                org_hisat = FavoriteSpecies.FAVORITE_SPECIES[org]['hisat2']
                org_gtf = FavoriteSpecies.FAVORITE_SPECIES[org]['gtf']
            except KeyError:
                DieGracefully.die_gracefully(
                    DieGracefully.BAD_ORG, a['organism'])
            a['hisat2_idx'] = org_hisat
            a['gtf'] = org_gtf
        # Convert all of the paths into absolute paths
        a['fq_folder'] = os.path.realpath(
            os.path.expanduser(str(a['fq_folder'])))
        a['hisat2_idx'] = os.path.realpath(
            os.path.expanduser(str(a['hisat2_idx'])))
        a['gtf'] = os.path.realpath(os.path.expanduser(str(a['gtf'])))
        a['outdir'] = os.path.realpath(os.path.expanduser(str(a['outdir'])))
        a['workdir'] = os.path.realpath(os.path.expanduser(str(a['workdir'])))
        if a['expr_groups']:
            a['expr_groups'] = os.path.realpath(os.path.expanduser(str(
                a['expr_groups'])))
        try:
            assert a['headcrop'] >= 0
            assert isinstance(a['headcrop'], int)
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--headcrop')
        try:
            assert a['mincts'] >= 0
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--min-cts')
        try:
            assert a['rrna_screen'] >= 0
            assert isinstance(a['rrna_screen'], int)
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--rrna_screen')
        try:
            assert a['subsample'] >= 0
            if a['subsample'] > 0:
                assert a['subsample'] >= a['rrna_screen']
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--subsample')
        try:
            assert a['mem'] >= 24000
            assert isinstance(a['mem'], int)
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--mem')
        try:
            assert a['tmp_space'] >= 0
            assert isinstance(a['tmp_space'], int)
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--tmp')
        try:
            assert a['walltime'] >= 2
            assert isinstance(a['walltime'], int)
        except AssertionError:
            DieGracefully.die_gracefully(
                DieGracefully.BAD_NUMBER, '--walltime')
        self.pipe_logger.debug('GTF: %s', a['gtf'])
        self.pipe_logger.debug('Adapters: %s', a['adapters'])
        self.pipe_logger.debug('FASTQ Folder: %s', a['fq_folder'])
        self.pipe_logger.debug('Output Dir: %s', a['outdir'])
        self.pipe_logger.debug('Working Dir: %s', a['workdir'])
        self.pipe_logger.debug('HISAT2 Idx: %s', a['hisat2_idx'])
        self.pipe_logger.debug('Expr Groups: %s', a['expr_groups'])
        self.pipe_logger.debug('Strandness: %s', a['strand'])
        # Check that the adapters and GTF file exist
        try:
            handle = open(a['gtf'], 'rt')
            handle.close()
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_GTF)
        if not a['adapters']:
            a['adapters'] = '/home/msistaff/public/CHURP_Deps/v1/db/all_illumina_adapters.fa'
        else:
            try:
                a['adapters'] = os.path.realpath(
                    os.path.expanduser(str(a['adapters'])))
                handle = open(a['adapters'], 'rt')
                handle.close()
            except OSError:
                DieGracefully.die_gracefully(DieGracefully.BAD_ADAPT)
        # Validate the FASTQ folder
        self._validate_fastq_folder(a['fq_folder'])
        # Validate the hisat2 index
        self._validate_hisat_idx(a['hisat2_idx'])
        # Validate the experimental groups xlsx 
        a = self._validate_groupsheet(a)
        # Sanitize the hisat2 index path
        a['hisat2_idx'] = dir_funcs.sanitize_path(
            a['hisat2_idx'], self.pipe_logger)
        return a

    def _validate_fastq_folder(self, d):
        """Raise an error if the FASTQ directory does not exist or does not
        have any FASTQ files."""
        try:
            contents = os.listdir(d)
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_FASTQ)
        # Check if there is at least one file ending in a standard fastq suffix
        fq_pat = re.compile(r'^.+((.fq(.gz)?$)|(.fastq(.gz)?$))')
        has_fastq = False
        for f in contents:
            if re.match(fq_pat, f):
                has_fastq = True
                break
        if has_fastq:
            return
        else:
            DieGracefully.die_gracefully(DieGracefully.EMPTY_FASTQ)
        return

    def _validate_hisat_idx(self, i):
        """Raise an error if the provided HISAT2 index is not complete -
        all of the [1-8].ht2l? files should be present."""
        # Build glob patterns for the normal and long indices
        norm_idx = i + '.[1-8].ht2'
        long_idx = i + '.[1-8].ht2l'
        # Do the search
        self.pipe_logger.debug('Searching for %s', norm_idx)
        norm_idx_files = glob.glob(norm_idx)
        self.pipe_logger.debug('Found %i idx files', len(norm_idx_files))
        # There should be 8 total
        if len(norm_idx_files) == 8:
            return
        else:
            self.pipe_logger.debug(
                'Normal idx not found. Searching for long idx.')
            long_idx_files = glob.glob(long_idx)
            self.pipe_logger.debug(
                'Found %i long idx files', len(long_idx_files))
            if len(long_idx_files) == 8:
                return
            else:
                self.pipe_logger.error('Cound not find HISAT2 idx files!')
                DieGracefully.die_gracefully(DieGracefully.BAD_HISAT)
        return

    def _prepare_samplesheet(self):
        """Call the samplesheet build method here. This will build the
        dictionary that will hold all samplesheet data, and then write it into
        the output directory."""
        self._run_checks()
        is_pe = self.sheet.compile(self.real_out, self.real_work)
        # We want to throw an error if there is a mix of PE and SE samples.
        if len(set(is_pe)) > 1:
            # If we get here, then we should separate the list of samples into
            # those that are SE and those that are PE
            se = []
            pe = []
            for samp in sorted(self.sheet.final_sheet):
                # Check for entries in the R1 and R2 slots
                r1 = self.sheet.final_sheet[samp]['FastqR1files']
                r2 = self.sheet.final_sheet[samp]['FastqR2file']
                grp = self.sheet.final_sheet[samp]['Group']
                sname = samp + ' (Group: ' + grp + ')'
                if r1 and not r2:
                    se.append(sname)
                elif r1 and r2:
                    pe.append(sname)
            DieGracefully.die_gracefully(DieGracefully.PE_SE_MIX, pe, se)
        ss_path = self.sheet.write_sheet(self.real_out, self.pipe_name, '|')
        return ss_path
        
    def _validate_groupsheet(self, args):
        # make a stub if one was not created.
        # @return the args with accurate args['expr_groups'] 
        if args['expr_groups'] is None:
            args['expr_groups'] = self._create_stub_groupsheet(args)
        # Now do basic validation
        groups = args['expr_groups']
        # throw an error if the expr_groups file is not an excel file 
        if not groups[-4:] == "xlsx" and not groups[-3:] == "xls":
            self.sheet_logger.error(
                'A filetype different than an excel spreadsheet was '
                'supplied for --expr_groups. You supplied a filetype '
                f'end in {groups[-4:]}. Please replace this with an excel '
                'spreadsheet in which the first sheet has two columns '
                'titled "SampleName" and "Group"')
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_NO_SAMP_GPS)                 
        # load the group sheet (the first sheet) from the xlsx
        groups_sheet = pd.read_excel(groups, 
                                     sheet_name = 0, 
                                     keep_default_na = False)
        # throw an error if SampleName is not a header in the group sheet
        if not 'SampleName' in groups_sheet:
            self.sheet_logger.error(
                'The xlsx experimental groups first sheet must have a '
                'column named "SampleName"')
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_NO_SAMP_GPS)   
        # throw an error if Group is not a header in the group sheet
        if not 'Group' in groups_sheet:
            self.sheet_logger.error(
                'The xlsx experimental groups first sheet must have a '
                'column named "Group"')
            DieGracefully.die_gracefully(DieGracefully.BRNASEQ_NO_SAMP_GPS)
        return(args)

    def _create_stub_groupsheet(self, args):
        ## This produces an experimental groups xlsx with SampleNames
        ## It requires us to figure out the sample names, which means I'm
        ## duplicating abit of code from SampleSheet
        samples = self._get_sample_names(args['fq_folder'])
        samples.sort()
        
        # populate the group sheet with the sample names
        groups = pd.DataFrame({'SampleName': samples,
                    'Group': 'NULL'})
        contrasts = pd.DataFrame({'Comparison_Name' : [],
                      'Reference_Group' : [],
                      'Test_Group' : []})
        # make sure the directory exists and make the expt group path
        expr_group_path = f'{args["outdir"]}/experimental_groups.xlsx'
        # if the outdir hasn't been made, make it
        if not os.path.isdir(args['outdir']):
            os.makedirs(args['outdir'])
        # save the group and contrast sheets
        with pd.ExcelWriter(expr_group_path) as writer:  
            groups.to_excel(writer, sheet_name='groups', index = False)
            contrasts.to_excel(writer, sheet_name='contrasts', index = False)
        # return the new path to add back to args
        return(os.path.realpath(expr_group_path))

    def _get_sample_names(self, fq_dir):
        # Four different regexs to determine if its a fastq then grab the 
        # sample name, for fastqs from UMGC vs. SRA
        samp_re = re.compile(
            r'(_S[0-9]+)?'
            r'(_[ATCG]{4,})?'
            r'(_L00[1-8])?'
            r'(_R(1|2))?_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        # RE to get forwards reads fq files
        fq_re = re.compile(
            r'^.+[^_R2]_001\.((fq(\.gz)?$)|(fastq(\.gz)?$))',
            flags=re.I)
        sra_re = re.compile(r'^.+_1\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        sra_samp_re = re.compile(r'_(1|2)\.((fq(\.gz)?$)|(fastq(\.gz)?$))')
        # iterate through contents of the fastq dir and grab the sample names
        fq_dir_contents = os.listdir(fq_dir)
        sample_names = []
        for current_file in fq_dir_contents:
            # if it is a fastq from UMGC...
            if re.match(fq_re, current_file):
                # Extract the samplename from the fastq name
                sn = re.sub(samp_re, '', current_file)
                sample_names.append(sn)
            # if it is a fastq from SRA...
            elif re.match(sra_re, current_file):
                sn = re.sub(sra_samp_re, '', current_file)
                sample_names.append(sn)
        return(sample_names)


    def qsub(self):
        """Write the qsub command. We will need the path to the samplesheet,
        the number of samples in the samplesheet, and the scheduler options
        that were passed as arguments. This is defined in the subclass rather
        than in the main Pipeline class because the exact form of the qsub
        command depends on which pipeline we are running."""
        gs = self.valid_args["expr_groups"]
        ss = self._prepare_samplesheet()
        # Make the qsub array key
        keyname = default_files.default_array_key(self.pipe_name)
        keyname = os.path.join(self.real_out, keyname)
        if os.path.isfile(keyname):
            self.pipe_logger.warning(
                'Sbatch key file %s exists. Overwriting!', keyname)
        try:
            handle = open(keyname, 'wt')
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        # The sheet is sorted in this way before it is written to disk, so it
        # should be safe to sort it this way
        handle.write('Sbatch.Index\tSampleName\n')
        for index, samplename in enumerate(sorted(self.sheet.final_sheet)):
            handle.write(str(index+1) + '\t' + samplename + '\n')
        handle.flush()
        handle.close()
        # Make the script filename
        pname = default_files.default_pipeline(self.pipe_name)
        pname = os.path.join(self.real_out, pname)
        if os.path.isfile(pname):
            self.pipe_logger.warning(
                'Submission script %s already exists. Overwriting!', pname)
        try:
            handle = open(pname, 'wt')
        except OSError:
            DieGracefully.die_gracefully(DieGracefully.BAD_OUTDIR)
        # Write the header of the script
        handle.write('#!/bin/bash\n')
        # Write some meta-data lines
        handle.write('# Generated by CHURP version ' +
                     CHURPipelines.__version__ + '\n')
        handle.write('# Generated at ' + CHURPipelines.NOW + '\n')
        handle.write('set -e\n')
        handle.write('set -u\n')
        handle.write('set -o pipefail\n')
        self.pipe_logger.debug(
            'Number of samples: %i', len(self.sheet.final_sheet))
        self.pipe_logger.debug('Samplesheet: %s', ss)
        # Write command to figure out email address of submitting user
        handle.write('user_name="$(id -u -n)"\n')
        handle.write('user_email="${user_name}@umn.edu"\n')
        # Set the group string here
        if self.group:
            qsub_group = '-A ' + self.group
        else:
            qsub_group = ''
        qsub_array = '1'
        if len(self.sheet.final_sheet) > 1:
            qsub_array += '-' + str(len(self.sheet.final_sheet))
        # Write a few variables into the header of the script so they are
        # easy to find
        handle.write('CHURP_VERSION=' + '"' + CHURPipelines.__version__ + '"\n')
        handle.write('SUMMARY_ONLY=' + '"' + self.summary_only + '"\n')
        handle.write('KEYFILE=' + '"' + keyname + '"\n')
        handle.write('QSUB_ARRAY=' + '"' + qsub_array + '"\n')
        handle.write('OUTDIR=' + '"' + str(self.real_out) + '"\n')
        handle.write('WORKDIR=' + '"' + str(self.real_work) + '"\n')
        handle.write('DE_SCRIPT=' + '"' + self.de_script + '"\n')
        handle.write('REPORT_SCRIPT=' + '"' + self.report_script + '"\n')
        handle.write('GROUPSHEET=' + '"' + gs + '"\n')
        handle.write('SAMPLESHEET=' + '"' + ss + '"\n')
        handle.write('PURGE=' + '"' + self.purge + '"\n')
        handle.write('RRNA_SCREEN=' + '"' + self.rrna_screen + '"\n')
        handle.write('SUBSAMPLE=' + '"' + self.subsample + '"\n')
        handle.write('PIPE_SCRIPT="$(cd "$( dirname "${BASH_SOURCE[0]}" )" '
                     '>/dev/null && pwd )/$(basename $0)"\n')
        # These are the variables we want to export into the single sample job
        # script.
        single_cmd_vars = ','.join([
            'SampleSheet="${SAMPLESHEET}"',
            'PURGE="${PURGE}"',
            'RRNA_SCREEN="${RRNA_SCREEN}"',
            'SUBSAMPLE="${SUBSAMPLE}"'
            ])
        # And the command for the single-sample job array
        aln_cmd = [
            'sbatch',
            '--parsable',
            '--ignore-pbs',
            '-p', self.msi_queue,
            '--mail-type=BEGIN,END,FAIL',
            '--mail-user="${user_email}"',
            qsub_group,
            '-o', '"${OUTDIR}/bulk_rnaseq_single_sample-%A.%a.out"',
            '-e', '"${OUTDIR}/bulk_rnaseq_single_sample-%A.%a.err"',
            '-N', '1',
            '--mem=' + str(self.mem) + 'mb',
            '--tmp=' + str(self.tmp_space) + 'mb',
            '-n', '1',
            '-c', str(self.ppn),
            '--time=' + str(self.walltime * 60),
            '--array="${QSUB_ARRAY}"',
            '--export=' + single_cmd_vars,
            self.single_sample_script,
            '||',
            'exit',
            '1']
        # The variables to export for the summary job
        summary_vars = ','.join([
            'SampleSheet="${SAMPLESHEET}"',
            'GroupSheet="${GROUPSHEET}"',
            'CHURP_VERSION="${CHURP_VERSION}"',
            'MINLEN="' + self.min_gene_len + '"',
            'MINCPM="' + self.min_cts + '"',
            'RSUMMARY="${DE_SCRIPT}"',
            'PIPE_SCRIPT="${PIPE_SCRIPT}"',
            'BULK_RNASEQ_REPORT="${REPORT_SCRIPT}"'])
        # Write some logic to detect if we are running in a job allocation.
        # Quit the bash script if it is.
        handle.write('if [ ! -z "${SLURM_JOB_ID+NULL}" ]\n')
        handle.write('    then echo "You should run this script with \'bash\' from outside of a job allocation." > /dev/stderr\n')
        handle.write('    exit 99\n')
        handle.write('fi\n')
        # Make two versions of the summary command: one with the single sample
        # array dependency and one without
        summary_cmd_dep = [
            'sbatch',
            '--parsable',
            '--ignore-pbs',
            '-p', self.msi_queue,
            '--mail-type=BEGIN,END,FAIL',
            '--mail-user="${user_email}"',
            qsub_group,
            '-o', '"${OUTDIR}/run_summary_stats-%j.out"',
            '-e', '"${OUTDIR}/run_summary_stats-%j.err"',
            '-N', '1',
            '--mem=' + str(self.mem) + 'mb',
            '--tmp=' + str(self.tmp_space) + 'mb',
            '-n', '1',
            '-c', str(self.ppn),
            '--time=' + str(self.walltime * 60),
            '--depend=afterok:${single_id}',
            '--export=' + summary_vars,
            self.summary_script,
            '||',
            'exit',
            '1']
        summary_cmd = [
            'sbatch',
            '--parsable',
            '--ignore-pbs',
            '-p', self.msi_queue,
            '--mail-type=BEGIN,END,FAIL',
            '--mail-user="${user_email}"',
            qsub_group,
            '-o', '"${OUTDIR}/run_summary_stats-%j.out"',
            '-e', '"${OUTDIR}/run_summary_stats-%j.err"',
            '-N', '1',
            '--mem=' + str(self.mem) + 'mb',
            '--tmp=' + str(self.tmp_space) + 'mb',
            '-n', '1',
            '-c', str(self.ppn),
            '--time=' + str(self.walltime * 60),
            '--export=' + summary_vars,
            self.summary_script,
            '||',
            'exit',
            '1']
        # Write the logic for if the user has specified the "--summary-only"
        # option.
        handle.write('if [ -f "${OUTDIR}/.in_progress" ]\n')
        handle.write('then\n')
        handle.write('    echo "Job already submitted; aborting." > /dev/stderr\n')
        handle.write('    exit 1\n')
        handle.write('fi\n')
        handle.write('if [ "${SUMMARY_ONLY}" = "true" ]\n')
        handle.write('then\n')
        handle.write('    touch "${OUTDIR}/.in_progress"\n')
        handle.write('    summary_id=$(' + ' '.join(summary_cmd) + ')\n')
        handle.write('else\n')
        handle.write('    touch "${OUTDIR}/.in_progress"\n')
        handle.write('    single_id=$(' + ' '.join(aln_cmd) + ')\n')
        handle.write('    summary_id=$(' + ' '.join(summary_cmd_dep) + ')\n')
        handle.write('fi\n')
        # Write some echo statements for users' information
        handle.write('echo "You are running CHURP version ${CHURP_VERSION}"\n')
        handle.write('echo "Output and logs will be written to ${OUTDIR}"\n')
        handle.write('echo "Emails will be sent to ${user_email}"\n')
        handle.write('echo "Sbatch array to samplename key: ${KEYFILE}"\n')
        handle.write('if [ "${SUMMARY_ONLY}" = "true" ]\n')
        handle.write('    then echo "--summary-only" specified. No single samples job array ID\n')
        handle.write('    else echo "Single samples job array ID: ${single_id}"\n')
        handle.write('fi\n')
        handle.write('echo "Summary job ID: ${summary_id}"\n')
        self.pipe_logger.debug('sbatch:\n%s', ' '.join(aln_cmd))
        self.pipe_logger.debug('sbatch:\n%s', ' '.join(summary_cmd))
        handle.flush()
        handle.close()
        # Check if we want to automatically submit the script
        if self.nosubmit:
            qsub_dat = None
        else:
            qsub_cmd = ['bash', pname]
            qsub_proc = subprocess.Popen(
                qsub_cmd,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            qsub_stdout, qsub_stderr = qsub_proc.communicate()
            qsub_dat = (qsub_stdout, qsub_stderr, qsub_proc)
        return (pname, ss, keyname, qsub_dat)
