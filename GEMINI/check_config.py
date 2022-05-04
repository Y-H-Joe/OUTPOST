#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 20:30:40 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####

=================================== input =====================================

=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================

####=======================================================================####
"""
import sys
from collections import Counter
# import re
import os

def check_config(df_config):
    """ check the config.tsv. """
    def rm_nan_from_set(s):
        return {x for x in s if x==x}

    # the column names should be formated
    samples_col,fq_dir_col,bam_dir_col,assembly_col,assembly_dir_col,group_col = 0,0,0,0,0,0
    cols = df_config.columns

    if len(cols) > len(rm_nan_from_set(set(cols))):
        sys.exit("GEMINI: detected duplicate columns in contig.tsv. exit.")

    for col in cols:
        if col == 'samples':
            samples_col += 1
            continue

        elif col == 'fq_dir':
            fq_dir_col += 1
            continue

        elif col == 'bam_dir':
            bam_dir_col += 1
            continue

        elif col == 'assembly':
            assembly_col += 1
            continue

        elif col == 'assembly_dir':
            assembly_dir_col += 1
            continue

        elif col == 'group':
            group_col += 1
            continue

    assert all([x == 1 for x in [samples_col,fq_dir_col,bam_dir_col,assembly_col,assembly_dir_col,group_col]]),\
    "GEMINI: detected errors in columns of contig.tsv. exit."

    # the id of each sample should be unique
    samples = df_config['samples']
    assert len(samples) == len(rm_nan_from_set(set(samples))),\
    "GEMINI: sample ids are not unique. exit."
    sample_list = [x.strip() for x in samples]

    # the bam of each sample should be unique
    bams = df_config['bam_dir']
    assert len(bams) == len(rm_nan_from_set(set(bams))),\
    "GEMINI: bam directory is not unique. exit."
    bam_list = [x.strip() for x in bams]

    # the fq of each sample should be unique
    fqs = df_config['fq_dir']
    assert len(fqs) == len(rm_nan_from_set(set(fqs))),\
    "GEMINI: fastq files are not unique. exit."
    fq_list = [x.strip() for x in fqs]

    # only one assembly is accepted each run
    assert len(rm_nan_from_set(set(df_config['assembly']))) == 1,\
    "GEMINI: multiple assemblies. exit."
    assembly = list(df_config['assembly'])[0]

    # same assembly should have the same dir
    assert len(rm_nan_from_set(set(df_config['assembly_dir']))) == 1,\
    "GEMINI: multiple assembly_dir. exit."
    assembly_dir = list(df_config['assembly_dir'])[0]

    # each group should have at least 3 samples
    # return comparison_dict:{'group1':[sample1,sample3],'group2':[sample2,sample3]}
    groups_list = df_config['group']
    assert len(set(groups_list)) == len(rm_nan_from_set(set(groups_list))),\
    "GEMINI: detected sample which belongs to no group. exit."
    groups_list = [x.strip().split(',') for x in groups_list]
    group_list = [val for sublist in groups_list for val in sublist]
    group_count = dict(Counter(group_list))
    assert all([x >= 3 for x in group_count.values()]),\
    "GEMINI: each group should have at least 3 samples. exit."
    comparison_dict = {}
    for group in group_count.keys():
        samples = []
        for index,groups_of_this_sample in enumerate(groups_list):
            if group in groups_of_this_sample:
                samples.append(sample_list[index])
        comparison_dict[group] = samples

    # detect files existence
    files = bam_list + fq_list + [assembly_dir]
    assert all([os.path.exists(x) for x in files]),\
    "GEMINI: detect non-exist files (fq/bam/assembly). exit."

    # detect unwanted symbols
    assert "," not in "".join(sample_list + fq_list + bam_list) + assembly_dir + assembly,\
    "GEMINI: ',' is not allowed in names. exit."

    # duplicate sample_list/fq_list/bam_list/comparison_dict for downstream comparison analysis
    if len(group_list) == 1:
        print("GEMINI: only one group detected. Duplicate to smooth downstream analysis")
        sample_list_copy = [sample + '_copy' for sample in sample_list]
        sample_list = sample_list + sample_list_copy
        fq_list = fq_list * 2
        bam_list = bam_list * 2
        comparison_dict.update({(group_list[0] + '_copy') : sample_list_copy})

    return sample_list, fq_list, bam_list, assembly, assembly_dir, comparison_dict


