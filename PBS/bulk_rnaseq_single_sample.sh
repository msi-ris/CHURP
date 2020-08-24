#!/bin/bash -l

set -e
set -u
set -o pipefail

# Export the PS4 variable for the trace
# Taken from https://wiki.bash-hackers.org/scripting/debuggingtips
LOG_SECTION="General"
export PS4='+[$(date "+%F %T")] [${SLURM_JOB_ID}] [${LOG_SECTION}]: '

# Define a function to report Slurm performance metrics to the log files.
# slurm_res_report() {
    # Use sstat to get the resource usage info
    # sacct -X -j "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}" -l
    # sstat -a -j "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    # SACCT_STATS=$(sacct -X "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}" -P -n --format=Elapsed)
    # SSTAT_STATS=$(sstat -j "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.batch" -P -n --format=MaxDiskRead,MaxDiskWrite,MaxRSS)
    # # Chew up the disk read/write values
    # DISK_READ_BYTES=$(echo "${SSTAT_STATS}" | cut -f 1 -d '|')
    # DISK_WRITE_BYTES=$(echo "${SSTAT_STATS}" | cut -f 2 -d '|')
    # DISK_READ_GB=$(echo "scale=3; ${DISK_READ_BYTES}/1024/1024/1024" | bc -l)
    # DISK_WRITE_GB=$(echo "scale=3; ${DISK_WRITE_BYTES}/1024/1024/1024" | bc -l)
    # # # And print them out:
    # echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Job resource summary"
    # echo "Job ID: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    # echo "MSI Group: ${SLURM_JOB_ACCOUNT}"
    # echo "Partition (queue): ${SLURM_JOB_PARTITION}"
    # echo "Number of CPUs: ${SLURM_NPROCS}"
    # echo "Approx. Walltime: $(echo ${SACCT_STATS} | cut -f 1 -d '|')"
    # echo "Memory Requested: ${SLURM_MEM_PER_NODE}"
    # echo "Approx. Memory Used: $(echo ${SSTAT_STATS} | cut -f 3 -d '|')"
    # echo "Approx. Data Read Off Disk: ${DISK_READ_GB}GB"
    # echo "Approx. Data Written to Disk: ${DISK_WRITE_GB}GB"
# }

# Define a function to report errors to the job log and give meawningful exit
# codes. This just wraps a bunch of exit calls into a case block
pipeline_error() {
    # Take the pipeline section as a positional argument
    case "${1}" in
    "General")
        echo "${SampleSheet} is incompatible with this version of CHURP." > /dev/stderr
        echo "${SampleSheet} was generated with version ${SAMPLESHEET_VERSION}, and this script requires ${PIPELINE_VERSION}." > /dev/stderr
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 100
        ;;
    "Subsampling")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "There was an error during read subsampling!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 101
        ;;
    "rRNA.Subsampling")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "There was an error during read subsampling for rRNA abundance estimation!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 102
        ;;
    "BBDuk")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "BBDuk encountered an error!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 103
        ;;
    "FastQC.Raw")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "FastQC failed on the raw FASTQ files!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 104
        ;;
    "Trimmomatic")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "Trimmomatic encountered an error!" >> "${LOG_FNAME}"
        echo "Please see the error message above for details." >> "${LOG_FNAME}"
        echo "If you see a message about 'Error: Unable to detect quality encoding' then your FASTQ files do not have a standard quality encoding." >> "${LOG_FNAME}"
        echo "If you see a Java exception and you have specified custom Trimmomatic options, then this suggests a problem with your option string." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 105
        ;;
    "FastQC.Trimmed")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "FastQC failed on the trimmed FASTQ files!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 106
        ;;
    "HISAT2")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "HISAT2 encountered an error!" >> "${LOG_FNAME}"
        echo "If you specified custom options for HISAT2, then this indicates a problem with your options string." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 107
        ;;
    "MarkDuplicates")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "Picard MarkDuplicates encountered an error!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 108
        ;;
    "BAM.Filtering")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "SAMTools encountered an error while filtering the HISAT2 alignment!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 109
        ;;
    "BAM.Coord.Sort")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "SAMTools encountered an error while coordinate-sorting the HISAT2 alignment!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 110
        ;;
    "BAM.Stats")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "SAMTools encountered an error while generating summary statistics for the HISAT2 alignment!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 111
        ;;
    "InsertSizeMetrics")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "Picard InsertSizeMetrics encountered an error!" >> "${LOG_FNAME}"
        echo "Please see the error messages above for details." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 112
        ;;
    "Alignment.Summary")
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "CHURP encountered an error while extracting alignment summaries!" >> "${LOG_FNAME}"
        echo "This means the alignment summary output file from HISAT2 has gone missing." >> "${LOG_FNAME}"
        echo "Please re-run CHURP with the --purge option to rerun the pipeline from the beginning." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 113
        ;;
    *)
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "CHURP encountered an undefined error!" >> "${LOG_FNAME}"
        echo "Please send the CHURP command, version, samplesheet, and pipeline.sh script to help@msi.umn.edu for debugging." >> "${LOG_FNAME}"
        # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 200
        ;;
    esac
}

# Check the major version of the pipeline. If they mismatch, then quit with an error
PIPELINE_VERSION="0"
SAMPLESHEET_VERSION=$(tail -n 1 "${SampleSheet}" | sed -E 's/#//g')
if [ "${SAMPLESHEET_VERSION}" -ne "${PIPELINE_VERSION}" ]
then
    pipeline_error "${LOG_SECTION}"
fi

# Parse SampleSheet using $SLURM_ARRAY_TASK_ID
# Assume there is a header line, so for ID n, the sample is n+1
# Sample Sheet: SampleNM, R1, R2, Trim (yes or no), TrimmomaticOption, Hisat2Index, Hisat2Option, GTF/GFF
IN=$(head -n "${SLURM_ARRAY_TASK_ID}" "${SampleSheet}" | tail -1)
# handle empty line and comment line. We return a 0 exit status because we do not
# want a comment/blank line in the samplesheet to hold up the array jobs
[ -z "${IN// }" ] && echo "You have submitted an empty line from the sample sheet to the job array. This array job will quit without error, but you should determine why this occured." && exit 0
[ "$IN" = "#.*" ] && echo "You have submitted a comment line from the sample sheet to the job array. This array job will quit without error, but you should determine why this occured." && exit 0
while IFS="|" read -ra OPTS; do 
    SAMPLENM=${OPTS[0]}
    EXPR_GROUP=${OPTS[1]}
    R1FILE=${OPTS[2]}
    R2FILE=${OPTS[3]}
    OUTDIR=${OPTS[4]}
    WORKDIR=${OPTS[5]}
    TRIM=${OPTS[6]}
    RMDUP=${OPTS[7]}
    TRIMOPTS=${OPTS[8]}
    HISAT2INDEX=${OPTS[9]}
    HISAT2OPTS=${OPTS[10]}
    STRAND=${OPTS[11]}
    GTFFILE=${OPTS[12]}
done <<< "$IN"

# Start the trace. In this case, we use file descriptor 5 to avoid clobbering
# any other fds that are in use
LOGDIR="${OUTDIR}/Logs"
mkdir -p "${LOGDIR}"
TRACE_FNAME="${LOGDIR}/${SAMPLENM}_Trace.log"
LOG_FNAME="${LOGDIR}/${SAMPLENM}_Analysis.log"
# Write the samplename to the .e PBS file
echo "# $(date '+%F %T') Slurm error file for ${SAMPLENM}" >> /dev/stderr
echo "# $(date '+%F %T'): For a human-readable log, see ${LOG_FNAME}" >> /dev/stderr
echo "# $(date '+%F %T'): For a debugging trace, see ${TRACE_FNAME}" >> /dev/stderr
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
# We don't want to clobber the old files; just append to them
echo "###############################################################################" >> "${LOG_FNAME}"
echo "# $(date '+%F %T'): Analysis started for ${SAMPLENM}" >> "${LOG_FNAME}"
echo "# $(date '+%F %T'): Job ID: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}" >> "${LOG_FNAME}"
# This is annoying to see, so we will exclude it from verbose outout
# load necessary modules
module load python3/3.6.3_anaconda5.0.1
module load hisat2/2.1.0
module load fastqc/0.11.7
module load trimmomatic/0.33
module load samtools/1.7
module load picard/2.9.0
module load R/3.5.0
module load java/openjdk-11.0.2

module list -t
echo '#BEGIN_MODULES' >> "${LOG_FNAME}"
module list -t 2>> "${LOG_FNAME}"
echo '#END_MODULES' >> "${LOG_FNAME}"

# For future debugging, print which java we are using
echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Using $(which java)" >> "${LOG_FNAME}"

# Enable trace debugging here
exec 5>> "${TRACE_FNAME}"
echo "##### BEGIN TRACE FOR ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} #####" >> "${TRACE_FNAME}"
export BASH_XTRACEFD=5
set -x

# Set paths to BBDuk, seqtk, and the SILVA databases
DEPS_DIR="/home/msistaff/public/CHURP_Deps/v${PIPELINE_VERSION}"
SILVA_REF="${DEPS_DIR}/db/SILVA_132_LSU_SSU_Ref_Dedup_Kmers_min100.fasta.gz"
BBDUK="${DEPS_DIR}/Supp/bbmap/bbduk.sh"
SEQTK="${DEPS_DIR}/Supp/seqtk-1.3/seqtk"
GTFTOGENEPRED="${DEPS_DIR}/Supp/UCSC/gtfToGenePred"
CONDA_ENV="${DEPS_DIR}/Conda/envs/RNASeQC_RefPrep"
COLLAPSE_GTF="${DEPS_DIR}/Supp/GTEx_Pipeline/collapse_annotation.py"
RNASEQC="${DEPS_DIR}/Supp/RNASeQC/rnaseqc.v2.3.4.linux"

# Check if we are running in paired or single end mode
if [ -z "${R2FILE}" ]
then
    PE="false"
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): No R2 file detected; running ${SAMPLENM} as single-end" >> "${LOG_FNAME}"
else
    PE="true"
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): R2 file detected; running ${SAMPLENM} as paired-end" >> "${LOG_FNAME}"
fi

# check whether to purge files or not. $PURGE will be parsed by command line
if [ "${PURGE}" = "true" ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): PURGE=true; deleting work directory for ${SAMPLENM} and re-running all analyses." >> "${LOG_FNAME}"
    rm -rf "${WORKDIR}/singlesamples/${SAMPLENM}"
fi

# set working directory
mkdir -p "${WORKDIR}/singlesamples/${SAMPLENM}" && cd "${WORKDIR}/singlesamples/${SAMPLENM}"

# start workflow with check point
if [ -f "${SAMPLENM}.done" ]; then
    echo "Found completed analysis, exit" >> "${LOG_FNAME}"
    # slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
    exit 0
fi

echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Subsampling"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ "${SUBSAMPLE}" -eq 0 ]; then
    echo "# $(date '+%F %T'): Not subsampling reads for sample ${SAMPLENM} for analysis" >> "${LOG_FNAME}"
else
    echo "# $(date '+%F %T'): Subsampling ${SAMPLENM} to ${SUBSAMPLE} fragments for rRNA quantification" >> "${LOG_FNAME}"
    ${SEQTK} sample -s123 -2 "${R1FILE}" "${SUBSAMPLE}" | gzip -c > "${WORKDIR}/singlesamples/${SAMPLENM}/Subsample_R1.fastq.gz" 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    R1FILE="${WORKDIR}/singlesamples/${SAMPLENM}/Subsample_R1.fastq.gz"
    if [ "${PE}" = "true" ]; then
        ${SEQTK} sample -s123 -2 "${R2FILE}" "${SUBSAMPLE}" | gzip -c > "${WORKDIR}/singlesamples/${SAMPLENM}/Subsample_R2.fastq.gz" 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
        R2FILE="${WORKDIR}/singlesamples/${SAMPLENM}/Subsample_R2.fastq.gz"
    fi
fi


echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="rRNA.Subsampling"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f subsamp.done ]; then
    # subsample the FASTQ and assay for rRNA contamination
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Subsampling reads to ${RRNA_SCREEN} fragments." >> "${LOG_FNAME}"
    ${SEQTK} sample -s123 -2 "${R1FILE}" "${RRNA_SCREEN}" > "${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_R1.fastq" 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    if [ "${PE}" = "true" ]; then
        ${SEQTK} sample -s123 -2 "${R2FILE}" "${RRNA_SCREEN}" > "${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_R2.fastq" 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    fi
    touch subsamp.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found subsampled reads" >> "${LOG_FNAME}"
fi

echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="BBDuk"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
# Check if the BBDuk analysis has been finished
if [ ! -f bbduk.done ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Using BBDuk to search for rRNA contamination in subsampled reads." >> "${LOG_FNAME}"
    if [ "${PE}" = "true" ]; then
        ${BBDUK} \
            in="${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_R1.fastq" \
            in2="${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_R2.fastq" \
            ref="${SILVA_REF}" \
            stats="${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_rRNA_Stats.txt" \
            k=25 \
            editdistance=1 \
            prealloc=t \
            threads="${SLURM_NPROCS}" \
            -Xmx10g \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    else
        ${BBDUK} \
            in="${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_R1.fastq" \
            ref="${SILVA_REF}" \
            stats="${WORKDIR}/singlesamples/${SAMPLENM}/BBDuk_rRNA_Stats.txt" \
            k=25 \
            editdistance=1 \
            prealloc=t \
            threads="${SLURM_NPROCS}" \
            -Xmx10g \
             2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    fi
    touch bbduk.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found complete BBDuk analysis." >> "${LOG_FNAME}"
fi

echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="FastQC.Raw"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f fastqc.done ]; then
    if [ "${PE}" = "true" ]
    then
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running fastqc on ${R1FILE} and ${R2FILE}." >> "${LOG_FNAME}"
        fastqc \
            -t 2 \
            --extract \
            --outdir="${WORKDIR}/singlesamples/${SAMPLENM}" \
            "${R1FILE}" \
            "${R2FILE}" \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
        && touch fastqc.done
    else
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running fastqc on ${R1FILE}." >> "${LOG_FNAME}"
        fastqc \
            --extract \
            --outdir="${WORKDIR}/singlesamples/${SAMPLENM}" \
            "${R1FILE}" \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
            && touch fastqc.done
    fi
fi

# The fastqc files are auto-extracted, so we search for the directories that
# end in '_fastqc', pull out the total number of reads, then throw away the
# directory
for fastqc_out_dir in $(find "${WORKDIR}/singlesamples/${SAMPLENM}" -maxdepth 1 -type d -name '*_fastqc')
do
    read_no=$(basename "${fastqc_out_dir}" | sed -nr 's/.*_R?(1|2)_(001_)?fastqc/\1/p')
    if [ -z "${read_no}" ]
    then
        read_no="1"
    fi
    out_pref="${SAMPLENM}_${read_no}.raw_readcount.txt"
    echo "${SAMPLENM} $(grep '^Total Sequences' ${fastqc_out_dir}/fastqc_data.txt | cut -f 2)" > "${out_pref}"
    awk '/>>Per base sequence quality/{flag=1; next} />>END_MODULE/{flag=0} flag' \
        "${fastqc_out_dir}/fastqc_data.txt" \
        | sed -e 's/ /./g' \
        > "${WORKDIR}/singlesamples/${SAMPLENM}/${SAMPLENM}_${read_no}.raw_quals.txt"
    rm -rf "${fastqc_out_dir}.processed"
    mv -f "${fastqc_out_dir}" "${fastqc_out_dir}.processed"
done

# To use -basein option, fastq file must be in format *_R1_*.fastq and *_R2_*.fastq
# trimmomatic can handle the -basein and -baseout options if there is only one
# read (single end).
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Trimmomatic"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ "${TRIM}" = "yes" ]; then
    if [ ! -f trimmomatic.done ]; then
        if [ "${PE}" = "true" ]
        then
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running trimmomatic on ${R1FILE} and ${R2FILE}." >> "${LOG_FNAME}"
            java -jar "${TRIMMOMATIC}"/trimmomatic.jar \
                PE \
                -threads "${SLURM_NPROCS}" \
                "${R1FILE}" "${R2FILE}" \
                "${SAMPLENM}_1P.fq.gz" "${SAMPLENM}_1U.fq.gz" "${SAMPLENM}_2P.fq.gz" "${SAMPLENM}_2U.fq.gz" \
                $(echo "${TRIMOPTS}" | envsubst) \
                2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
                && touch trimmomatic.done
        else
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running trimmomatic on ${R1FILE}." >> "${LOG_FNAME}"
            java -jar "${TRIMMOMATIC}"/trimmomatic.jar \
                SE \
                -threads "${SLURM_NPROCS}" \
                "${R1FILE}" \
                "${SAMPLENM}_trimmed.fq.gz" \
                $(echo "${TRIMOPTS}" | envsubst) \
                2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
                && touch trimmomatic.done
        fi
    fi
    echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
    LOG_SECTION="FastQC.Trimmed"
    echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
    if [ ! -f fastqc.trim.done ]; then
        if [ "${PE}" = "true" ]
        then
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running fastqc on trimmed fastq files." >> "${LOG_FNAME}"
            fastqc \
                -t 2 \
                --extract \
                --outdir="${WORKDIR}/singlesamples/${SAMPLENM}" \
                "${SAMPLENM}_1P.fq.gz" \
                "${SAMPLENM}_2P.fq.gz" \
                2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
                && touch fastqc.trim.done
        else
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Running fastqc on trimmed fastq file." >> "${LOG_FNAME}"
            fastqc \
                --extract \
                --outdir="${WORKDIR}/singlesamples/${SAMPLENM}" \
                "${SAMPLENM}_trimmed.fq.gz" \
                2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}" \
                && touch fastqc.trim.done
        fi
    fi
fi

for fastqc_out_dir in $(find "${WORKDIR}/singlesamples/${SAMPLENM}" -maxdepth 1 -type d -name '*_fastqc')
do
    read_no=$(basename "${fastqc_out_dir}" | sed -nr 's/.*(1|2)P_fastqc/\1/p')
    # If the files are single-end, then $read_no is empty
    if [ -z "${read_no}" ]
    then
        read_no="1"
    fi
    out_pref="${SAMPLENM}_${read_no}.trimmed_readcount.txt"
    echo "${SAMPLENM} $(grep '^Total Sequences' ${fastqc_out_dir}/fastqc_data.txt | cut -f 2)" > "${out_pref}"
    awk '/>>Per base sequence quality/{flag=1; next} />>END_MODULE/{flag=0} flag' \
        "${fastqc_out_dir}/fastqc_data.txt" \
        | sed -e 's/ /./g' \
        > "${WORKDIR}/singlesamples/${SAMPLENM}/${SAMPLENM}_${read_no}.trim_quals.txt"
    rm -rf "${fastqc_out_dir}.processed"
    mv -f "${fastqc_out_dir}" "${fastqc_out_dir}.processed"
done

# HISAT2 chokes on quoted reads. We have to do this dumb quoting strategy because
# some filenames may have spaces in them, and this protects it. My thought is that
# the Perl wrapper script splits arguments with spaces in them
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="HISAT2"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f hisat2.done ]; then
    if [ "${TRIM}" = "yes" ]; then
        if [ "${PE}" = "true" ]
        then
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Aligning trimmed reads with HISAT2." >> "${LOG_FNAME}"
            # This string is so ugly because we we have to quote the arguments to hisat2 in a weird way to protect them from splitting
            hisat2 \
                ${HISAT2OPTS} \
                -x "${HISAT2INDEX}" \
                -1 <(gzip -cd "${SAMPLENM}_1P.fq.gz" || cat "${SAMPLENM}_1P.fq") \
                -2 <(gzip -cd "${SAMPLENM}_2P.fq.gz" || cat "${SAMPLENM}_2P.fq") \
                2> alignment.summary \
                | samtools view -hb -o "${SAMPLENM}.bam" - \
                && touch hisat2.done \
                || pipeline_error "${LOG_SECTION}"
        else
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Aligning trimmed reads with HISAT2." >> "${LOG_FNAME}"
            hisat2 \
                ${HISAT2OPTS} \
                -x "${HISAT2INDEX}" \
                -U <(gzip -cd "${SAMPLENM}_trimmed.fq.gz" || cat "${SAMPLENM}_trimmed.fq.gz") \
                2> alignment.summary \
                | samtools view -hb -o "${SAMPLENM}.bam" - \
                && touch hisat2.done \
                || pipeline_error "${LOG_SECTION}"
        fi
    else
        if [ "${PE}" = "true" ]
        then
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Aligning reads with HISAT2." >> "${LOG_FNAME}"
            hisat2 \
                ${HISAT2OPTS} \
                -x "${HISAT2INDEX}" \
                -1 <(gzip -cd "${R1FILE}" || cat "${R1FILE}") \
                -2 <(gzip -cd "${R2FILE}" || cat "${R2FILE}") \
                2> alignment.summary \
                | samtools view -hb -o "${SAMPLENM}.bam" - \
                && touch hisat2.done \
                || pipeline_error "${LOG_SECTION}"
        else
            echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Aligning reads with HISAT2." >> "${LOG_FNAME}"
            hisat2 \
                ${HISAT2OPTS} \
                -x "${HISAT2INDEX}" \
                -U <(gzip -cd "${R1FILE}" || cat "${R1FILE}") \
                2> alignment.summary \
                | samtools view -hb -o "${SAMPLENM}.bam" - \
                && touch hisat2.done \
                || pipeline_error "${LOG_SECTION}"
        fi
    fi
    # Stick the alignment summary onto the analysis log
    cat alignment.summary >> "${LOG_FNAME}"
fi

# Next, mark or remove duplicates
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="MarkDuplicates"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f dup.done ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Soring raw HISAT2 BAM by query in prep for deduplication." >> "${LOG_FNAME}"
    _JAVA_OPTIONS="-Djava.io.tmpdir=${WORKDIR}/singlesamples/${SAMPLENM}/picard_tmp" ${PTOOL}/picard.jar \
        SortSam \
        I="${SAMPLENM}.bam" \
        O="${SAMPLENM}_Raw_QuerySort.bam" \
        SORT_ORDER="queryname" \
        2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    if [ "${RMDUP}" = "yes" ]; then 
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Removing duplicate reads with Picard MarkDuplicates." >> "${LOG_FNAME}"
        _JAVA_OPTIONS="-Djava.io.tmpdir=${WORKDIR}/singlesamples/${SAMPLENM}/picard_tmp" ${PTOOL}/picard.jar \
            MarkDuplicates \
            I="${SAMPLENM}_Raw_QuerySort.bam" \
            O="${SAMPLENM}_DeDup.bam" \
            REMOVE_DUPLICATES="true" \
            ASSUME_SORT_ORDER="queryname" \
            M="${SAMPLENM}_MarkDup_Metrics.txt" \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
        # Set the name of the bam to filter
        TO_FLT="${SAMPLENM}_DeDup.bam"
    else
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Marking duplicate reads with Picard MarkDuplicates." >> "${LOG_FNAME}"
        _JAVA_OPTIONS="-Djava.io.tmpdir=${WORKDIR}/singlesamples/${SAMPLENM}/picard_tmp" ${PTOOL}/picard.jar \
            MarkDuplicates \
            I="${SAMPLENM}_Raw_QuerySort.bam" \
            O="${SAMPLENM}_MarkDup.bam" \
            REMOVE_DUPLICATES="false" \
            ASSUME_SORT_ORDER="queryname" \
            M="${SAMPLENM}_MarkDup_Metrics.txt" \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
        TO_FLT="${SAMPLENM}_MarkDup.bam"
    fi
    touch dup.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found deduplicated/marked BAM files." >> "${LOG_FNAME}"
    # But, be sure to set the TO_FLT variable:
    if [ "${RMDUP}" = "yes" ]; then
        TO_FLT="${SAMPLENM}_DeDup.bam"
    else
        TO_FLT="${SAMPLENM}_MarkDup.bam"
    fi
fi

# Next, remove unmapped reads and reads with MAPQ<60
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="BAM.Filtering"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f mapq_flt.done ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Removing unmapped and MAPQ<60 reads for counting." >> "${LOG_FNAME}"
    samtools view \
        -bhu \
        -@ "${SLURM_NPROCS}" \
        -F 4 \
        -q 60 \
        -o "${SAMPLENM}_Filtered.bam" \
        "${TO_FLT}"
    touch mapq_flt.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found filtered BAM for counting." >> "${LOG_FNAME}"
fi

# Next, sort by coord for IGV purposes
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="BAM.Coord.Sort"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f coord_sort.done ]; then
    echo "$ ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T' ): Removing any old SAMtools sort files." >> "${LOG_FNAME}"
    find "${WORKDIR}/singlesamples/${SAMPLENM}" \
        -mindepth 1 \
        -maxdepth 1 \
        -regextype posix-extended \
        -regex '.*/temp\.[0-9]{4}+.bam' \
        -exec rm {} \;
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Sorting filtered BAM file by coordinate." >> "${LOG_FNAME}"
    samtools sort \
        -O bam \
        -@ "${SLURM_NPROCS}" \
        -T temp \
        -o "${SAMPLENM}_Filtered_CoordSort.bam" \
        "${SAMPLENM}_Filtered.bam" \
        2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Sorting raw BAM file by coordinate." >> "${LOG_FNAME}"
    samtools sort \
        -O bam \
        -@ "${SLURM_NPROCS}" \
        -T temp \
        -o "${SAMPLENM}_Raw_CoordSort.bam" \
        "${TO_FLT}" \
        2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Indexing coordinate-sorted BAM files." >> "${LOG_FNAME}"
    samtools index "${SAMPLENM}_Filtered_CoordSort.bam"
    samtools index "${SAMPLENM}_Raw_CoordSort.bam"
    touch coord_sort.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found sorted and indexed BAM files." >> "${LOG_FNAME}"
fi

# Set these variables for linking at the end of the script
FOR_COUNTS="${SAMPLENM}_Filtered.bam"
RAW_COORD="${SAMPLENM}_Raw_CoordSort.bam"
RAW_COORD_IDX="${SAMPLENM}_Raw_CoordSort.bam.bai"
FLT_COORD="${SAMPLENM}_Filtered_CoordSort.bam"
FLT_COORD_IDX="${SAMPLENM}_Filtered_CoordSort.bam.bai"

# Generate some stats on the raw BAM for the report
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="BAM.Stats"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f bamstats.done ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Generating alignment stats based on raw BAM." >> "${LOG_FNAME}"
    samtools stats "${RAW_COORD}" > "${SAMPLENM}_bamstats.txt"
    touch bamstats.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found raw BAM stats." >> "${LOG_FNAME}"
fi

# Try the RNASeQC metrics gathering
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="RNASeQC"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
echo "# $(date '+%F %T'): Note, this section is OPTIONAL (errors will not kill pipeline jobs)." >> /dev/stderr
if [ ! -f rnaseqc.done ]; then
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): 'Collapsing' gene models in GTF for use with RNASeQC." >> "${LOG_FNAME}"
    source activate "${CONDA_ENV}" || true
    python "${COLLAPSE_GTF}" <(gzip -cd "${GTFFILE}" || cat "${GTFFILE}") "${WORKDIR}/singlesamples/${SAMPLENM}/collapsed.gtf" || true
    source deactivate || true
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Collecting unstranded RNAseq metrics with RNASeQC." >> "${LOG_FNAME}"
    RNASEQC_OPTIONS="-v -v --sample=${SAMPLENM}_Unstranded --legacy"
    "${RNASEQC}" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/collapsed.gtf" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/${RAW_COORD}" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/RNASeQC_Out" \
        ${RNASEQC_OPTIONS} 2>> "${LOG_FNAME}" || true
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Collecting stranded RNAseq metrics with RNASeQC." >> "${LOG_FNAME}"
    RNASEQC_OPTIONS="-v -v --sample=${SAMPLENM} --legacy"
    # A bit strange - if the data are single-read data, then the strand has to
    # be flipped. This was confirmed with single-read pico v2 and pico v1 data.
    if [[ "${STRAND}" = "2" && "${PE}" = "true" ]]; then
        RNASEQC_OPTIONS="${RNASEQC_OPTIONS} --stranded=RF"
    elif [[ "${STRAND}" = "1" && "${PE}" = "true" ]]; then
        RNASEQC_OPTIONS="${RNASEQC_OPTIONS} --stranded=FR"
    elif [[ "${STRAND}" = "2" && "${PE}" = "false" ]]; then
        RNASEQC_OPTIONS="${RNASEQC_OPTIONS} --stranded=FR"
    elif [[ "${STRAND}" = "1" && "${PE}" = "false" ]]; then
        RNASEQC_OPTIONS="${RNASEQC_OPTIONS} --stranded=RF"
    fi
    "${RNASEQC}" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/collapsed.gtf" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/${RAW_COORD}" \
        "${WORKDIR}/singlesamples/${SAMPLENM}/RNASeQC_Out" \
        ${RNASEQC_OPTIONS} 2>> "${LOG_FNAME}" || true
    touch rnaseqc.done
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found RNAseq metrics checkpoint." >> "${LOG_FNAME}"
fi

# Use picard to collect the insert size metrics, but only if paired end
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="InsertSizeMetrics"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ ! -f is_stats.done ]; then
    if [ "${PE}" = "true" ]; then
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Collecting insert size metrics with Picard InsertSizeMetrics." >> "${LOG_FNAME}"
        mkdir -p "${OUTDIR}/InsertSizeMetrics"
        mkdir -p "${WORKDIR}/singlesamples/${SAMPLENM}/picard_tmp"
        _JAVA_OPTIONS="-Djava.io.tmpdir=${WORKDIR}/singlesamples/${SAMPLENM}/picard_tmp" ${PTOOL}/picard.jar \
            CollectInsertSizeMetrics \
            I="${WORKDIR}/singlesamples/${SAMPLENM}/${FOR_COUNTS}" \
            O="${OUTDIR}/InsertSizeMetrics/${SAMPLENM}_metrics.txt" \
            H="${OUTDIR}/InsertSizeMetrics/${SAMPLENM}_hist.pdf" \
            2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}"
        # Extract the mean, median, standard deviation from the metrics file
        if [ -s "${OUTDIR}/InsertSizeMetrics/${SAMPLENM}_metrics.txt" ]; then
            grep \
                -A 1 \
                '^MEDIAN_INSERT_SIZE' \
                "${OUTDIR}/InsertSizeMetrics/${SAMPLENM}_metrics.txt" \
                | tail -n 1 \
                | cut -f 1,5,6,9,11,15,17 \
                > "${WORKDIR}/singlesamples/${SAMPLENM}/IS_Stats.txt"
        else
            # Echo seven NA into the IS stats file
            echo -e 'NA\tNA\tNA\tNA\tNA\tNA\tNA' \
            > "${WORKDIR}/singlesamples/${SAMPLENM}/IS_Stats.txt"
        fi
    else
        echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Sample is single-read. No insert size metrics possible." >> "${LOG_FNAME}"
    fi
else
    echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Found insert size metrics." >> "${LOG_FNAME}"
fi

# Use awk to pick apart the alignment summary
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Alignment.Summary"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
if [ "${PE}" = "true" ]
then
    TOTAL_READS=$(awk '/Total pairs:/ {print $NF}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    UNMAP=$(awk '/Aligned concordantly or discordantly 0 time/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    SINGLE_MAP=$(awk '/Aligned concordantly 1 time/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    MULTI_MAP=$(awk '/Aligned concordantly >1 times/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    DISCO_MAP=$(awk '/Aligned discordantly 1 time/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
else
    TOTAL_READS=$(awk '/Total reads:/ {print $NF}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    UNMAP=$(awk '/Aligned 0 time/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    SINGLE_MAP=$(awk '/Aligned 1 time/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    MULTI_MAP=$(awk '/Aligned >1 times/ {F=NF-1; print $F}' alignment.summary 2>> "${LOG_FNAME}" || pipeline_error "${LOG_SECTION}")
    DISCO_MAP="NA"
fi
echo "${SAMPLENM} ${TOTAL_READS} ${UNMAP} ${SINGLE_MAP} ${MULTI_MAP} ${DISCO_MAP}" > "hisat_map_summary.txt"

# the final step is to link the sorted.rmdup.bam file to the allsamples/
# directory, as just the samplename. This is a bit of a hack to get featureCounts
# to not print huge paths as samplenames
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Cleanup"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Linking ${SAMPLENM} into ${WORKDIR}/allsamples" >> "${LOG_FNAME}"
mkdir -p "${WORKDIR}/allsamples"
ln -sf "${WORKDIR}/singlesamples/${SAMPLENM}/${FOR_COUNTS}" "${WORKDIR}/allsamples/${SAMPLENM}"

echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Linking coordinate-sorted BAMs into ${OUTDIR}/Coordinate_Sorted_BAMs/" >> "${LOG_FNAME}"
mkdir -p "${OUTDIR}/Coordinate_Sorted_BAMs"
ln -sf "${WORKDIR}/singlesamples/${SAMPLENM}/${RAW_COORD}" "${OUTDIR}/Coordinate_Sorted_BAMs/${RAW_COORD}"
ln -sf "${WORKDIR}/singlesamples/${SAMPLENM}/${RAW_COORD_IDX}" "${OUTDIR}/Coordinate_Sorted_BAMs/${RAW_COORD_IDX}"
ln -sf "${WORKDIR}/singlesamples/${SAMPLENM}/${FLT_COORD}" "${OUTDIR}/Coordinate_Sorted_BAMs/${FLT_COORD}"
ln -sf "${WORKDIR}/singlesamples/${SAMPLENM}/${FLT_COORD_IDX}" "${OUTDIR}/Coordinate_Sorted_BAMs/${FLT_COORD_IDX}"

echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Finished processing ${SAMPLENM}." >> "${LOG_FNAME}"
touch "${SAMPLENM}.done"

# Finally, let's clean up
echo "# ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} $(date '+%F %T'): Removing HISAT2 bam, markdup/dedup bam, and raw querysort bam to reduce disk usage." >> "${LOG_FNAME}"
rm -f "${SAMPLENM}.bam" "${SAMPLENM}_MarkDup.bam" "${SAMPLENM}_DeDup.bam" "${SAMPLENM}_Raw_QuerySort.bam"

# And, run our stats function
# slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr

# And close the file descriptor we were using for the trace
exec 5>&-
