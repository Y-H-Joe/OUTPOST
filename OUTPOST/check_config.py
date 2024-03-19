#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 20:30:40 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
df_config: dataframe
rm_batch_effect: true or false
# mode: reads | assembly | bam
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
    samples_col,refgenome_col, r1_fq_dir_col, r2_fq_dir_col, se_fq_dir_col,bam_dir_col,\
        assembly_col,assembly_dir_col,group_col,batch_col = 0,0,0,0,0,0,0,0,0,0
    cols = df_config.columns

    if len(cols) > len(rm_nan_from_set(set(cols))):
        sys.exit("OUTPOST: detected duplicate columns in contig.tsv. exit.")

    for col in cols:
        if col == 'samples':
            samples_col += 1
        elif col == 'refgenome':
            refgenome_col += 1
        elif col == 'r1_fq_dir':
            r1_fq_dir_col += 1
        elif col == 'r2_fq_dir':
            r2_fq_dir_col += 1
        elif col == 'se_fq_dir':
            se_fq_dir_col += 1
        elif col == 'bam_dir':
            bam_dir_col += 1
        elif col == 'assembly':
            assembly_col += 1
        elif col == 'assembly_dir':
            assembly_dir_col += 1
        elif col == 'group':
            group_col += 1
        elif col == 'batch':
            batch_col += 1
    
    assert all([x == 1 for x in [samples_col,refgenome_col, r1_fq_dir_col, 
                                 r2_fq_dir_col, se_fq_dir_col, bam_dir_col, 
                                 assembly_col, assembly_dir_col, group_col, batch_col]]),\
    f"""OUTPOST: detected errors in columns names of {[samples_col,refgenome_col, r1_fq_dir_col, 
                                 r2_fq_dir_col, se_fq_dir_col, bam_dir_col, 
                                 assembly_col, assembly_dir_col, group_col, batch_col]}. exit."""

    ###================= check each column =================###
    # id
    # the id of each sample should be unique 
    samples = df_config['samples']
    assert len(samples) == len(rm_nan_from_set(set(samples))),\
    "OUTPOST: sample ids are not unique. exit."
    sample_list = [x.strip() for x in samples]

    # refgenome
    # the refgenome of each sample should match samples
    refgenome = df_config['refgenome']
    if len(rm_nan_from_set(refgenome)) > 0:
        assert len(refgenome) == len(samples),\
        "OUTPOST: the refgenome records not match sample records. exit."
        refgenome_list = [x.strip() for x in refgenome]
        for file in refgenome_list:
            if not os.path.exists(file):
                sys.exit(f"OUTPOST: detect {file} doesn't exist. exit")
    else:
        refgenome_list = []
        
    # fq
    # 保证fastq全有或者全不有。全有时，保证每一个sample都有。每个sample都有时，保证PE全有或全缺，若缺，则有SE顶上
    # 允许PE、SE混杂；允许不同sample有相同的fq
    # the fq of each sample should be unique
    r1_fqs = df_config['r1_fq_dir']
    r2_fqs = df_config['r2_fq_dir']
    se_fqs = df_config['se_fq_dir']
    # samples should all have fq or all not have fq
    if len(rm_nan_from_set(r1_fqs)) > 0 or len(rm_nan_from_set(r2_fqs)) > 0:
        r1_fq_list, r2_fq_list = [],[]
        for x in r1_fqs:
            try:
                r1_fq_list.append(x.strip())
            except:
                r1_fq_list.append('')
        for x in r2_fqs:
            try:
                r2_fq_list.append(x.strip())
            except:
                r2_fq_list.append('')
    else:
        r1_fq_list = []
        r2_fq_list = []
    if len(rm_nan_from_set(se_fqs)) > 0:
        se_fq_list = []
        for x in se_fqs:
            try:
                se_fq_list.append(x.strip())
            except:
                se_fq_list.append('')
    else:
        se_fq_list = []

    if r1_fq_list.count('') == len(r1_fq_list):
        r1_fq_list = []
    if r2_fq_list.count('') == len(r2_fq_list):
        r2_fq_list = []
    if se_fq_list.count('') == len(se_fq_list):
        se_fq_list = []
        
    if len(r1_fq_list + r2_fq_list + se_fq_list) > 0: #samples should all have fq
        for file in r1_fq_list + r2_fq_list + se_fq_list:
            if file : # exclude ""
                if not os.path.exists(file):
                    sys.exit(f"OUTPOST: detect {file} doesn't exist. exit")
                
        for i,s in enumerate(sample_list):
            r1 = r1_fqs[i]
            r2 = r2_fqs[i]
            se = se_fqs[i]
            # r1 and r2 both exist or both nan
            assert (r1 == r1) is (r2 == r2), f"OUTPOST: sample {s}'s r1_fq, r2_fq doesn't match. exit."
            if not (r1 == r1): # if PE reads are nan, there must exist SE read
                assert (se == se), f"OUTPOST: sample {s} has no fastq. exit."
    
    
    # bam
    # 全有或全无；若全有，当存在
    # 允许不同sample有相同bam
    bams = df_config['bam_dir']
    bam_dir_list = []
    if len(rm_nan_from_set(bams)) > 0:
        for x in bams:
            try:
                bam_dir_list.append(x.strip())
            except:
                bam_dir_list.append('')
    if bam_dir_list.count('') == len(bam_dir_list):
        bam_dir_list = []
    if len(bam_dir_list) > 0:
        for file in bam_dir_list:
            if file : # exclude ""
                if not os.path.exists(file):
                    sys.exit(f"OUTPOST: detect {file} doesn't exist. exit")
    
    
    # assembly, assembly_dir
    # 需要与assembly_dir 同时有、同时没有
    # 允许部分有，但应该全部存在
    assemblys = df_config['assembly']
    assembly_dirs = df_config['assembly_dir']
    assembly_list, assembly_dir_list = [], []
    if len(rm_nan_from_set(assemblys)) > 0:
        for x in assemblys:
            try:
                assembly_list.append(x.strip())
            except:
                assembly_list.append('')
    if len(rm_nan_from_set(assembly_dirs)) > 0:
        for x in assembly_dirs:
            try:
                assembly_dir_list.append(x.strip())
            except:
                assembly_dir_list.append('')
    if assembly_list.count('') == len(assembly_list):
        assembly_list = []
    if assembly_dir_list.count('') == len(assembly_dir_list):
        assembly_dir_list = []
    for i,s in enumerate(sample_list):
        assembly = assemblys[i]
        assembly_dir = assembly_dirs[i]
        # assembly and assembly_dir both exist or both nan
        assert (assembly == assembly) is (assembly_dir == assembly_dir), f"OUTPOST: sample {s}'s assembly, assembly_dir doesn't match. exit."
    if len(assembly_dir_list) > 0: #samples should all have fq
        for file in assembly_dir_list:
            if file : # exclude ""
                if not os.path.exists(file):
                    sys.exit(f"OUTPOST: detect {file} doesn't exist. exit")


    # group
    # each group should have at least 3 samples
    # return comparison_dict:{'group1':[sample1,sample3],'group2':[sample2,sample3]}
    groups_list = df_config['group']
    assert len(set(groups_list)) == len(rm_nan_from_set(set(groups_list))),\
    "OUTPOST: detected sample which belongs to no group. exit."
    groups_list = [x.strip().split(',') for x in groups_list]
    group_list = [val.strip() for sublist in groups_list for val in sublist]
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
    
    # batch
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

    # detect unwanted symbols
    assert "," not in "".join(sample_list + refgenome_list + r1_fq_list + 
                              r2_fq_list + se_fq_list + bam_dir_list + assembly_list + assembly_dir_list), \
    "OUTPOST: ',' is not allowed in names. exit."

    return sample_list, refgenome_list, r1_fq_list, r2_fq_list, se_fq_list, bam_dir_list, assembly_list, assembly_dir_list, batch_list, comparison_dict

