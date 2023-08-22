#!/usr/bin/env python
"""Contains a dictionary for "favorite" species to be used as a shorthand for
specifying genome resources in the MSI bioref resource."""

import CHURPipelines

# The latest Ensembl releases pulled into bioref. Update this when bioref
# updtaes!
ENSEMBL_RELEASES = {
    'main': '109',
    'plants': '56',
    'fungi': '56',
    'metazoa': '56',
    'protists': '56'}

# This is potentially a bit silly, but it hopefully makes maintaining this
# easier. Store the ensembl division, species name (as encoded in the directory
# name) and the latest annotation/assembly version in a dictionary of tuples.
# These will have to be updated when bioref updates!!
# To add a new "favorite species" add an entry to this dictionary.
FAVE_ASM = {
    'human': ('main', 'Homo_sapiens', 'GRCh38.p13'),
    'mouse': ('main', 'Mus_musculus', 'GRCm39'),
    'rat': ('main', 'Rattus_norvegicus', 'mRatBN7.2'),
    'zebrafish': ('main', 'Danio_rerio', 'GRCz11'),
    'fly': ('main', 'Drosophila_melanogaster', 'BDGP6.32'),
    'worm': ('main', 'Caenorhabditis_elegans', 'WBcel235'),
    'yeast': ('main', 'Saccharomyces_cerevisiae', 'R64-1-1'),
    'cow': ('main', 'Bos_taurus', 'ARS-UCD1.2'),
    'dog': ('main', 'Canis_lupus_familiaris', 'ROS_Cfam_1.0'),
    'pig': ('main', 'Sus_scrofa', 'Sscrofa11.1'),
    'chicken': ('main', 'Gallus_gallus', 'bGalGal1.mat.broiler.GRCg7b'),
    'corn': ('plants', 'Zea_mays', 'Zm-B73-REFERENCE-NAM-5.0'),
    'soybean': ('plants', 'Glycine_max', 'Glycine_max_v2.1'),
    'barley': ('plants', 'Hordeum_vulgare',
        'MorexV3_pseudomolecules_assembly'),
    'wheat': ('plants', 'Triticum_aestivum', 'IWGSC'),
    'arabidopsis': ('plants', 'Arabidopsis_thaliana', 'TAIR10'),
    'chlamydomonas': ('plants', 'Chlamydomonas_reinhardtii',
        'Chlamydomonas_reinhardtii_v5.5'),
    'slime_mold': ('protists', 'Dictyostelium_discoideum', 'dicty_2.7')
    }

# Define the big data structure here.
FAVORITE_SPECIES = {}
for sp in FAVE_ASM:
    # Start an empty dictionary to hold the species data
    sp_d = {}
    # Unpack the tuple from our FAVE_ASM dictionary: get the ensembl division,
    # the species name, and the annotation version
    div, spn, ver = FAVE_ASM[sp]
    # Check the division - which Ensembl release should we use?
    rel = ENSEMBL_RELEASES.get(div)
    # Build the paths. Thankfully the names are programmatically generated, so
    # we can do this!!
    sp_d['bowtie'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/bowtie/genome')
    sp_d['bowtie2'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/bowtie2/genome')
    sp_d['bwa'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/bwa/genome')
    sp_d['hisat2'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/hisat2/genome')
    sp_d['seq'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/seq/genome.fa')
    sp_d['gtf'] = (
        f'{CHURPipelines.BIOREF_BASE}/{div}/{spn}-{rel}/{ver}/annotation/'
        f'{spn}.{ver}.{rel}.gtf')
    # Push this into the big dictionary
    FAVORITE_SPECIES[sp] = sp_d
