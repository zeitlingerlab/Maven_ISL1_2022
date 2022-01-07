"""
Author: Melanie Weilert
Affiliation: Stowers Institute
Aim: Pipeline for processing ChIP-seq data provided from the Srivastava lab.
Date: July 2021
Run: snakemake

Main target rules:
------------------

- bam_to_bw: Convert BAM reads to BW coverage files
    - seq: (1) deduplicate BAM file (2) resize fragments based on extension length (3) export to bw
- normalized_bw: Generate normalized .bw files
    - seq_bam_to_log2_normalized_seq_bw: Convert ChIP-seq BAM reads to normalized coverage files.
        - (1) deduplicate BAM file, (2) resize fragments to desired extension length, (3) normalize over WCE
- combined_bam_to_combined_norm_bw: Generate normalized .bw files across pooled sample
    - seq_bam_to_log2_normalized_seq_bw: Convert ChIP-seq BAM reads to normalized coverage files.
        - (1) deduplicate BAM file, (2) resize fragments to desired extension length, (3) normalize over WCE

"""

##########################################################################################
#Setup
##########################################################################################
import csv
import os
import math
import glob
import numpy as np
import pandas as pd
from itertools import product

##########################################################################################
#Import external data and prepare configuration parameters.
##########################################################################################
SAMPLES = pd.read_csv(f'tsv/samples.csv').set_index("sample_name", drop=False) #match genome with right indexes
SAMPLES.extension_length = SAMPLES.extension_length.fillna(0) #fill column with zeros for str(int(ext_length)) in macs2 rule
SAMPLES_WITH_NO_WCE = SAMPLES[~SAMPLES['sample_name'].str.contains('input')] #filter all rows that have "wce" in their name

##########################################################################################
# Request output files based on rules below
##########################################################################################
rule all:
    input:
        expand("data/bam/individual/{sample}.bam", sample = SAMPLES.sample_name),
        expand("data/bw/individual/{sample}.bw", zip, sample = SAMPLES.sample_name),
        expand("data/bam/combined/{sample}.bam", sample = SAMPLES.combined_sample_name.unique()),
        expand("data/bw/combined/{sample}.bw", zip, sample = SAMPLES.combined_sample_name.unique()),
        # expand("data/peaks/individual/{sample}_peaks.narrowPeak", zip, sample = SAMPLES_WITH_NO_WCE.output_sample),
        # expand("data/peaks/combined/{sample}_peaks.narrowPeak", zip, sample = SAMPLES_WITH_NO_WCE.combined_sample_name.unique()),
        # expand("data/idr/{sample}.log", zip, sample = SAMPLES_WITH_NO_WCE.sample_name),
        #
        # #Run snakemake -j 1 twice to trigger this last rule, dependent on idr files existing in the first place. :P
        # expand("data/idr/{idr_run}.bed", idr_run = [os.path.basename(s).replace('.txt', '') for s in glob.glob('data/idr/*.txt')]),

##########################################################################################
# Combination rules
##########################################################################################

rule bam_to_combined_bam:
    input:
        lambda wildcards: ['data/bam/individual/' + s + '.bam' for s in SAMPLES.sample_name[SAMPLES.combined_sample_name == wildcards.sample]]
    output:
        "data/bam/combined/{sample}.bam"
    message: "Merging .bam files:"
    shell:
        "samtools merge {output} {input}"

rule bam_to_norm_bw:
    input:
        wce = lambda wildcards: 'data/bam/individual/' + str(SAMPLES[SAMPLES.sample_name==wildcards.sample].control_name[0]) + '.bam',
        exp = "data/bam/individual/{sample}.bam",
    output:
        raw = "data/bw/individual/{sample}.bw",
        norm = "data/bw/individual/{sample}_log2_norm.bw"
    params:
        extension_length = lambda wildcards: str(int(SAMPLES[SAMPLES.sample_name==wildcards.sample].extension_length[0])),
        threads = 4,
    message: "Computing coverage of .bam file..."
    shell:
        """
        samtools index {input.wce}
        samtools index {input.exp}

        bamCoverage -b {input.exp} -o {output.raw} -of bigwig -bs 20 -p {params.threads} \
        -e {params.extension_length} --ignoreDuplicates --centerReads

        bamCompare -b1 {input.exp} -b2 {input.wce} -o {output.norm} -of bigwig \
        -e {params.extension_length} --ignoreDuplicates --centerReads \
        --scaleFactorsMethod readCount --operation log2 -bs 20 -p {params.threads}
        """

rule combined_bam_to_combined_norm_bw:
    input:
        exp = "data/bam/combined/{sample}.bam",
        wce = lambda wildcards: 'data/bam/combined/' + str(SAMPLES[SAMPLES.combined_sample_name==wildcards.sample].combined_control_name[0]) + '.bam',
    output:
        raw = "data/bw/combined/{sample}.bw",
        norm = "data/bw/combined/{sample}_log2_norm.bw"
    params:
        extension_length = lambda wildcards: str(int(SAMPLES[SAMPLES.combined_sample_name==wildcards.sample].extension_length[0])),
        threads = 4,
    message: "Computing coverage of .bam file..."
    shell:
        """
        samtools index {input.wce}
        samtools index {input.exp}

        bamCoverage -b {input.exp} -o {output.raw} -of bigwig -bs 20 -p {params.threads} \
        -e {params.extension_length} --ignoreDuplicates --centerReads

        bamCompare --bamfile1 {input.exp} --bamfile2 {input.wce} -o {output.norm} -of bigwig \
        -e {params.extension_length} --ignoreDuplicates --centerReads \
        --scaleFactorsMethod readCount --operation log2 -bs 20 -p {params.threads}
        """
