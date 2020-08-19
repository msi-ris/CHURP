#!/bin/bash

set -e
set -u
set -o pipefail

# Export the PS4 variable for the trace
# Taken from https://wiki.bash-hackers.org/scripting/debuggingtips
LOG_SECTION="General"
export PS4='+[$(date "+%F %T")] [${SLURM_JOB_ID}] [${LOG_SECTION}]: '

slurm_res_report() {
    # Use sstat to get the max VM size (for memory)
    MEM=$(sstat -j "${SLURM_JOB_ID}" -P -n --format=MaxVMSize)
    # Use sacct to get the other values:
    #   - Job name
    #   - Account
    #   - Used CPUTime
    #   - Used Walltime
    #   - Requested memory
    STATS=$(sacct -j "${SLURM_JOB_ID}.0" --format=JobName,Account,CPUTime,Elapsed,ReqMem -P -n)
    # And print them out:
    echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Job summary"
    echo "Job ID: ${SLURM_JOB_ID}"
    echo "Job Name: $(echo ${STATS} | cut -f 1 -d '|')"
    echo "MSI Group: $(echo ${STATS} | cut -f 2 -d '|')"
    echo "Number of CPUs: ${SLURM_NPROCS}"
    echo "CPU Time: $(echo ${STATS} | cut -f 3 -d '|')"
    echo "Walltime: $(echo ${STATS} | cut -f 4 -d '|')"
    echo "Memory Requested: $(echo ${STATS} | cut -f 5 -d '|')"
    echo "Memory Used: ${MEM}"
}

# Define a function to report errors to the job log and give meawningful exit
# codes. This just wraps a bunch of exit calls into a case block
pipeline_error() {
    # Take the pipeline section as a positional argument
    case "${1}" in
    "General")
        echo "${SampleSheet} is incompatible with this version of CHURP." > /dev/stderr
        echo "${SampleSheet} was generated with version ${SAMPLESHEET_VERSION}, and this script requires ${PIPELINE_VERSION}." > /dev/stderr
        slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 100
        ;;
    "featureCounts")
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "The featureCounts program failed to produce a counts matrix for your dataset." >> "${LOG_FNAME}"
        echo "Please check that you have permission and sufficient space to write to the output directory." >> "${LOG_FNAME}"
        slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 114
        ;;
    "edgeR")
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "The edgeR gene expresssion analysis script encountered an error." >> "${LOG_FNAME}"
        echo "Please check that your sample sheet is properly formatted and that your group labels do not contain special characters." >> "${LOG_FNAME}"
        echo "You may also find additonal R error messages below:" >> "${LOG_FNAME}"
        echo "" >> "${LOG_FNAME}"
        cat Rout.txt >> "${LOG_FNAME}"
        echo "" >> "${LOG_FNAME}"
        slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 115
        ;;
    "HTML.Report")
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "CHURP was unable to produce a summary HTML report for your run." >> "${LOG_FNAME}"
        echo "Your counts files and alignments should still be available." >> "${LOG_FNAME}"
        echo "Please send your pipeline.sh and samplesheet to help@msi.umn.edu for assistance." >> "${LOG_FNAME}"
        slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 119
        ;;
    *)
        echo "" >> "${LOG_FNAME}"
        echo "#### CHURP caught an error #####" >> "${LOG_FNAME}"
        echo "CHURP encountered an undefined error!" >> "${LOG_FNAME}"
        echo "Please send the CHURP command, version, samplesheet, and pipeline.sh script to help@msi.umn.edu for debugging." >> "${LOG_FNAME}"
        slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr
        exit 200
        ;;
    esac
}

# Check for PBS/Samplesheet version agreement
PIPELINE_VERSION="0"
SAMPLESHEET_VERSION=$(tail -n 1 "${SampleSheet}" | sed -E 's/#//g')
if [ "${SAMPLESHEET_VERSION}" -ne "${PIPELINE_VERSION}" ]
then
    pipeline_error "${LOG_SECTION}"
fi

# Read the first line of the samplesheet to get the working directory and the
# GTF. This is not super elegant, but it works with the samplesheet format
IN=$(head -n 1 "${SampleSheet}")
# handle empty line and comment line. We return a 0 exit status because we do not
# want a comment/blank line in the samplesheet to hold up the array jobs
[ -z "${IN// }" ] && echo "hit an empty line in sample sheet" && exit 0
[ "$IN" = "#.*" ] && echo "hit a comment line in sample sheet" && exit 0
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

# load necessary modules
module load R/3.5.0
module load python/3.6.3
# Make directories for the output files. We want to keep the outdir organized
LOGDIR="${OUTDIR}/Logs"
COUNTSDIR="${OUTDIR}/Counts"
PLOTSDIR="${OUTDIR}/Plots"
DEGDIR="${OUTDIR}/DEGs"
LOG_FNAME="${LOGDIR}/BulkRNASeq_Analysis.log"
TRACE_FNAME="${LOGDIR}/BulkRNASeq_Trace.log"

mkdir -p "${LOGDIR}" "${COUNTSDIR}" "${PLOTSDIR}" "${DEGDIR}"

# Write the samplename to the .e file
echo "# $(date '+%F %T') Standard error file for summary job" >> /dev/stderr
echo "# $(date '+%F %T'): For a human-readable log, see ${LOG_FNAME}" >> /dev/stderr
echo "# $(date '+%F %T'): For a debugging trace, see ${TRACE_FNAME}" >> /dev/stderr
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
# Print the loaded modules to stderr and a log file
module list -t
# We don't want to clobber the old files; just append to them
echo "###############################################################################" >> "${LOG_FNAME}"
echo "# $(date '+%F %T'): Summary job started" >> "${LOG_FNAME}"
echo "# $(date '+%F %T'): Job ID: ${SLURM_JOB_ID}" >> "${LOG_FNAME}"
echo '#BEGIN_MODULES' >> "${LOG_FNAME}"
module list -t 2>> "${LOG_FNAME}"
echo '#END_MODULES' >> "${LOG_FNAME}"

# Set up trace logging after loading modules to avoid dumping tons of module-related messages to the log
exec 5>> "${TRACE_FNAME}"
echo "##### BEGIN TRACE FOR ${SLURM_JOB_ID} #####" >> "${TRACE_FNAME}"
export BASH_XTRACEFD=5
set -x

if [ -z "${R2FILE}" ]
then
    PE="false"
    echo "# No R2 file detected; running as single-end" >> "${LOG_FNAME}"
    echo "# NOTE! We are assuming that all samples were subject to the same preparation and sequencing protocols." >> "${LOG_FNAME}"
    echo "# NOTE! We do not support mixing of single-end and paired-end samples for counts and differential expression testing." >> "${LOG_FNAME}"
else
    PE="true"
    echo "# R2 file detected; running as paired-end" >> "${LOG_FNAME}"
    echo "# NOTE! We are assuming that all samples were subject to the same preparation and sequencing protocols." >> "${LOG_FNAME}"
    echo "# NOTE! We do not support mixing of single-end and paired-end samples for counts and differential expression testing." >> "${LOG_FNAME}"
fi

# Set the R_LIBS_USER variable here. This is where packages will be loaded from
# within R
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Adding /home/msistaff/public/CHURP_Deps/v${PIPELINE_VERSION}/R to R_LIBS_USER" >> "${LOG_FNAME}"
export R_LIBS_USER="/home/msistaff/public/CHURP_Deps/v${PIPELINE_VERSION}/R"

# Set the path to the featureCounts executable.
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Using counting reads with featureCounts v. 1.6.2" >> "${LOG_FNAME}"
FEATURECOUNTS="/home/msistaff/public/CHURP_Deps/v${PIPELINE_VERSION}/Supp/subread-1.6.2-Linux-x86_64/bin/featureCounts"
mkdir -p "${WORKDIR}/allsamples" && cd "${WORKDIR}/allsamples"

# Use featureCounts to make a merged counts matrix
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="featureCounts"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Making a counts matrix for all samples." >> "${LOG_FNAME}"
BAM_LIST=($(find . -type l -exec basename {} \;| sort -V))
if [ "${PE}" = "true" ]
then
    echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Library is paired-end with strand ${STRAND}." >> "${LOG_FNAME}"
    "${FEATURECOUNTS}" \
        -a "${GTFFILE}" \
        -T ${SLURM_NPROCS} \
        -B \
        -p \
        -Q 10 \
        -s "${STRAND}" \
        -o subread_counts.txt \
        "${BAM_LIST[@]}" || pipeline_error "${LOG_SECTION}"
else
    echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Library is single-end with strand ${STRAND}" >> "${LOG_FNAME}"
    "${FEATURECOUNTS}" \
        -a "${GTFFILE}" \
        -T ${SLURM_NPROCS} \
        -Q 10 \
        -s "${STRAND}" \
        -o subread_counts.txt \
        "${BAM_LIST[@]}"   || pipeline_error "${LOG_SECTION}"
fi

# Summarize the merged count data, including descriptive summaries and differential expression tests if >1 group present.
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="edgeR"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Running edgeR analysis on counts" >> "${LOG_FNAME}"
Rscript \
    "${RSUMMARY}" \
    "${OUTDIR}" \
    "${WORKDIR}" \
    "${SampleSheet}" \
    "${WORKDIR}/allsamples/subread_counts.txt" \
    "${MINLEN}" \
    "${MINCPM}" \
    &> Rout.txt || pipeline_error "${LOG_SECTION}"

echo "# ----- Output from ${RSUMMARY} below" >> "${LOG_FNAME}"
cat Rout.txt >> "${LOG_FNAME}"
echo "# ----- End output from ${RSUMMARY}" >> "${LOG_FNAME}"

# Copy the merged counts into the output directory. The -u option to cp causes
# the copy to happen only if the source file is newer than the destination file
# or if the destination does not exist.
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="GTF.Summary"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Copying merged counts matrix and summary into ${COUNTSDIR}" >> "${LOG_FNAME}"
cp -u subread_counts.txt "${COUNTSDIR}/subread_counts.txt"
cp -u subread_counts.txt.summary "${COUNTSDIR}/subread_counts.txt.summary"

# We also want to keep the sorted BAM files and the GTF used for counts, in
# case the user wants to go back to it
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Copying GTF into ${OUTDIR}." >> "${LOG_FNAME}"
cp -u "${GTFFILE}" "${OUTDIR}"

# AWK command from S. Munro to make a translation table fo Ensembl IDs and
# gene names
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Generating translation table of gene name and Ensembl ID." >> "${LOG_FNAME}"
# awk -F '\t' '$3 == "gene" { print $9 }' <(gzip -cd "${GTFFILE}" || cat "${GTFFILE}") | tr -d ';"' | awk -F ' ' -v OFS='\t' '{print $2,$6}' > "${OUTDIR}/gene_id_gene_name_map.txt"
# This is a slower but more general command to get gene_id and gene_name
awk -F '\t' '$3 == "gene" {print $9}' <(gzip -cd "${GTFFILE}" || cat "${GTFFILE}") \
    | awk -F ';' -v OFS='\t' '{
    found_gene="false"
    for(i = 1; i <= NF; i++) {
        if( $i ~ /gene_id/ ) {
            split($i,GI," ")
        }
        if( $i ~ /gene_name/ ) {
            split($i,GN," ")
            found_gene="true"
        }
    }
    if(found_gene=="true") {
        print GI[2],GN[2]
    } else {
        print GI[2],GI[2]
    }
}' \
    | sed -e 's/"//g' \
    > "${OUTDIR}/gene_id_gene_name_map.txt"

# Chop up the subread_counts.txt file a little bit to make it easier to import
# into other tools like CLC Genomics Workbench
echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Linking gene symbol with Ensembl ID for subread_counts.txt" >> "${LOG_FNAME}"
echo -e "GeneName\t$(head -n 2 ${COUNTSDIR}/subread_counts.txt | tail -n 1 | cut -f 7-)" > "${COUNTSDIR}/subread_counts_gene_symbol.txt"
join \
    -e "NA" \
    -a 2 \
    <(sort -k 1,1 "${OUTDIR}/gene_id_gene_name_map.txt") \
    <(tail -n +3 "${COUNTSDIR}/subread_counts.txt" | cut -f 1,7- | sort -k 1,1) \
    | cut -f 2- -d ' ' \
    | tr ' ' '\t' \
    >> "${COUNTSDIR}/subread_counts_gene_symbol.txt"

# Link the work directories to the output directory
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Linking"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
rm -f "${OUTDIR}/singlesamples_work_directory" "${OUTDIR}/allsamples_work_directory"
ln -sf "${WORKDIR}/singlesamples" "${OUTDIR}/singlesamples_work_directory"
ln -sf "${WORKDIR}/allsamples" "${OUTDIR}/allsamples_work_directory"

# Prepare to generate HTML report by making unified files for everything. We
# will have to iterate over every sample
# Remove old summary files
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="Collate.Samples"
rm -f \
    "${WORKDIR}/allsamples/Read_Counts.txt" \
    "${WORKDIR}/allsamples/Samtools_Stats.txt" \
    "${WORKDIR}/allsamples/HISAT_Stats.txt" \
    "${WORKDIR}/allsamples/IS_Stats.txt" \
    "${WORKDIR}/allsamples/RNASeq_Metrics.txt" \
    "${WORKDIR}/allsamples/RNASeq_Metrics_Unstranded.txt" \
    "${WORKDIR}/allsamples/rRNA_kmers.txt"
# Use find to dig up all singlesample directories in the work director
SDIRS=($(find "${WORKDIR}/singlesamples" -maxdepth 1 -mindepth 1 -type d | sort -V))
for sampledir in "${SDIRS[@]}"
do
    sample=$(basename "${sampledir}")
    if [ ! -f "${sampledir}/${sample}.done" ]; then
        echo "# $(date '+%F %T'): No .done file for ${sample}, skipping." >> "${LOG_FNAME}"
        continue
    fi
    echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Collecting summaries for ${sample}" >> "${LOG_FNAME}"
    # Make the read counts files
    RAW_R1=$(find "${WORKDIR}/singlesamples/${sample}" -maxdepth 1 -type f -name '*1.raw_readcount.txt')
    RAW_R2=$(find "${WORKDIR}/singlesamples/${sample}" -maxdepth 1 -type f -name '*2.raw_readcount.txt')
    R1_COUNT=$(awk '{print $NF}' "${RAW_R1}")
    if [ -z "${RAW_R2}" ]
    then
        R2_COUNT="NA"
    else
        R2_COUNT=$(awk '{print $NF}' "${RAW_R2}")
    fi
    TRIMMED_R1=$(find "${WORKDIR}/singlesamples/${sample}" -maxdepth 1 -type f -name '*1.trimmed_readcount.txt')
    TRIMMED_R2=$(find "${WORKDIR}/singlesamples/${sample}" -maxdepth 1 -type f -name '*2.trimmed_readcount.txt')
    if [ -z "${TRIMMED_R1}" ]
    then
        TR1_COUNT="NA"
    else
        TR1_COUNT=$(awk '{print $NF}' "${TRIMMED_R1}")
    fi
    if [ -z "${TRIMMED_R2}" ]
    then
        TR2_COUNT="NA"
    else
        TR2_COUNT=$(awk '{print $NF}' "${TRIMMED_R2}")
    fi
    echo -e "${sample}\t${R1_COUNT}\t${R2_COUNT}\t${TR1_COUNT}\t${TR2_COUNT}" >> "${WORKDIR}/allsamples/Read_Counts.txt"
    cat "${WORKDIR}/singlesamples/${sample}/hisat_map_summary.txt" >> "${WORKDIR}/allsamples/HISAT_Stats.txt"
    # Then, get alignment summary files
    BAMSTAT="${WORKDIR}/singlesamples/${sample}/${sample}_bamstats.txt"
    NMAPPED=$(grep '^SN' "${BAMSTAT}" | grep 'raw total sequences' | cut -f 3)
    DUPLICATED=$(grep '^SN' "${BAMSTAT}" | grep 'reads duplicated' | cut -f 3)
    MQ0=$(grep '^SN' "${BAMSTAT}" | grep 'reads MQ0' | cut -f 3)
    MAXLEN=$(grep '^SN' "${BAMSTAT}" | grep 'maximum length' | cut -f 3)
    AVGLEN=$(grep '^SN' "${BAMSTAT}" | grep 'average length' | cut -f 3)
    AVGQUAL=$(grep '^SN' "${BAMSTAT}" | grep 'average quality' | cut -f 3)
    echo -e "${sample}\t${NMAPPED}\t${DUPLICATED}\t${MQ0}\t${MAXLEN}\t${AVGLEN}\t${AVGQUAL}" >> "${WORKDIR}/allsamples/Samtools_Stats.txt"
    # Get the BBDuk rRNA screen stats
    echo -e "${sample}\t$(grep '#Matched' ${WORKDIR}/singlesamples/${sample}/BBDuk_rRNA_Stats.txt | cut -f 2-3)" >> "${WORKDIR}/allsamples/rRNA_kmers.txt"
    # Then, get the insert size summaries. These are really simple.
    if [ -f "${WORKDIR}/singlesamples/${sample}/IS_Stats.txt" ]
    then
        echo -e "${sample}\t$(cat ${WORKDIR}/singlesamples/${sample}/IS_Stats.txt)" >> "${WORKDIR}/allsamples/IS_Stats.txt"
    fi
    # And collect the RNAseq metrics into a single file for summarization
    # Add the || true to the end of the command because these are optional steps
    sed -e "s/^/${sample}\t/g" <(tail -n +2 "${WORKDIR}/singlesamples/${sample}/RNASeQC_Out/${sample}.metrics.tsv") >> "${WORKDIR}/allsamples/RNASeq_Metrics.txt" || true
    sed -e "s/^/${sample}\t/g" <(tail -n +2 "${WORKDIR}/singlesamples/${sample}/RNASeQC_Out/${sample}_Unstranded.metrics.tsv") >> "${WORKDIR}/allsamples/RNASeq_Metrics_Unstranded.txt" || true
done

# We will generate an HTML report
# This is an ugly workaround, but sidesteps the problem of write-locked public dirs for Rmarkdown/knitr
echo "# $(date '+%F %T'): Finished section ${LOG_SECTION}" >> /dev/stderr
LOG_SECTION="HTML.Report"
echo "# $(date '+%F %T'): Entering section ${LOG_SECTION}" >> /dev/stderr
cp -u "${BULK_RNASEQ_REPORT}" "./Report.Rmd"
PATH=${PATH}:/panfs/roc/groups/14/msistaff/public/CHURP_Deps/v0/Supp/pandoc-2.3.1/bin Rscript -e "library(rmarkdown); rmarkdown::render('./Report.Rmd', output_file='"${OUTDIR}/Bulk_RNAseq_Report.html"', params=list(outdir='"${OUTDIR}"', workdir='"${WORKDIR}"', pipeline='"${PIPE_SCRIPT}"', samplesheet='"${SampleSheet}"'))" || pipeline_error "${LOG_SECTION}"

echo "# ${SLURM_JOB_ID} $(date '+%F %T'): Done summarizing bulk RNAseq run" >> "${LOG_FNAME}"
slurm_res_report | tee -a "${LOG_FNAME}" /dev/stderr

# Close the trace file descriptor
exec 5>&-
