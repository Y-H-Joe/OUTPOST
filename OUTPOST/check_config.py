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
import pandas as pd

def check_config(df_config, rm_batch_effect):
    if rm_batch_effect == 'None': rm_batch_effect = False

    """ check the config.tsv. """
    def rm_nan_from_set(s):
        return {x for x in s if x==x}

    # the column names should be formated
    samples_col,fq_dir_col,bam_dir_col,assembly_col,assembly_dir_col,group_col,batch_col = 0,0,0,0,0,0,0
    cols = df_config.columns

    if len(cols) > len(rm_nan_from_set(set(cols))):
        sys.exit("OUTPOST: detected duplicate columns in contig.tsv. exit.")

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

        elif col == 'batch':
            batch_col += 1
            continue

    assert all([x == 1 for x in [samples_col,fq_dir_col,bam_dir_col,assembly_col,assembly_dir_col,group_col,batch_col]]),\
    "OUTPOST: detected errors in columns of contig.tsv. exit."

    ###================= check each column =================###
    # the id of each sample should be unique
    samples = df_config['samples']
    assert len(samples) == len(rm_nan_from_set(set(samples))),\
    "OUTPOST: sample ids are not unique. exit."
    sample_list = [x.strip() for x in samples]

    # the bam of each sample should be unique
    bams = df_config['bam_dir']
    assert len(bams) == len(rm_nan_from_set(set(bams))),\
    "OUTPOST: bam directory is not unique. exit."
    bam_list = [x.strip() for x in bams]

    # the fq of each sample should be unique
    fqs = df_config['fq_dir']
    assert len(fqs) == len(rm_nan_from_set(set(fqs))),\
    "OUTPOST: fastq files are not unique. exit."
    fq_list = [x.strip() for x in fqs]

    # only one assembly is accepted each run
    assert len(rm_nan_from_set(set(df_config['assembly']))) == 1,\
    "OUTPOST: multiple assemblies. exit."
    assembly = list(df_config['assembly'])[0]

    # same assembly should have the same dir
    assert len(rm_nan_from_set(set(df_config['assembly_dir']))) == 1,\
    "OUTPOST: multiple assembly_dir. exit."
    assembly_dir = list(df_config['assembly_dir'])[0]

    # each group should have at least 3 samples
    # return comparison_dict:{'group1':[sample1,sample3],'group2':[sample2,sample3]}
    groups_list = df_config['group']
    assert len(set(groups_list)) == len(rm_nan_from_set(set(groups_list))),\
    "OUTPOST: detected sample which belongs to no group. exit."
    groups_list = [x.strip().split(',') for x in groups_list]
    group_list = [val for sublist in groups_list for val in sublist]
    group_count = dict(Counter(group_list))
    assert all([x >= 3 for x in group_count.values()]),\
    "OUTPOST: each group should have at least 3 samples. exit."
    comparison_dict = {}
    for group in group_count.keys():
        samples = []
        for index,groups_of_this_sample in enumerate(groups_list):
            if group in groups_of_this_sample:
                samples.append(sample_list[index])
        comparison_dict[group] = samples

    # batch column should be full or empty
    batch_list = [x for x in df_config['batch']]
    assert sum(pd.isna(batch_list)) in [0,len(batch_list)],\
        "OUTPOST: the batch column should be full or empty. exit."

    if sum(pd.isna(batch_list)) == 0:
        assert bool(rm_batch_effect) is bool(len(set(batch_list)) > 1),\
         "OUTPOST: chose not to 'rm_batch_effect' but detect multiple batches in the 'batch' column. Or chose to 'rm_batch_effect' but detect only 1 batch in the 'batch' column exit."
    if sum(pd.isna(batch_list)) == len(batch_list):
        assert bool(rm_batch_effect) is False,\
         "OUTPOST: No batch information in the 'batch' column, so only accept to not to 'rm_batch_effect'. exit."


    ###================= check others =================###
    # detect files existence
    files = bam_list + fq_list + [assembly_dir]
    for file in files:
        if not os.path.exists(file):
            sys.exit(f"OUTPOST: detect missing {file}. exit")
    # assert all([os.path.exists(x) for x in files]),\
    # "OUTPOST: detect missing files (fq/bam/assembly). exit."

    # detect unwanted symbols
    assert "," not in "".join(sample_list + fq_list + bam_list) + assembly_dir + assembly,\
    "OUTPOST: ',' is not allowed in names. exit."

    # duplicate sample_list/fq_list/bam_list/comparison_dict for downstream comparison analysis
    if len(group_list) == 1:
        print("OUTPOST: only one group detected. Duplicate to smooth downstream analysis")
        sample_list_copy = [sample + '_copy' for sample in sample_list]
        sample_list = sample_list + sample_list_copy
        fq_list = fq_list * 2
        bam_list = bam_list * 2
        comparison_dict.update({(group_list[0] + '_copy') : sample_list_copy})

    return sample_list, fq_list, bam_list, assembly, assembly_dir, comparison_dict


