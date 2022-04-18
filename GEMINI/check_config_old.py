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
import re
import os

def check_config(df_config):
    """ check the config.tsv. """
    def rm_nan_from_set(s):
        return {x for x in s if x==x}

    ## the column names should be formated
    samples_col,bam_dir_col,assembly_col,assembly_dir_col,group_col=\
    0,0,0,0,0

    assembly_list,assembly_loc_list, group_list, group_name_list = [],[],[],[]
    assembly_col_list,group_col_list=[],[]


    cols=df_config.columns

    if len(cols)>len(rm_nan_from_set(set(cols))):
        sys.exit("detected duplicate columns in contig.tsv. exit.")

    for col in cols:
        if col=='samples':
            samples_col+=1
            continue

        elif col=='bam_dir':
            bam_dir_col+=1
            continue

        elif re.match('assembly\d+$',col):
            assembly_col+=1
            assembly_dir_col+=1
            assembly_dir=col+"_dir"
            if not assembly_dir in cols:
                sys.exit("detected {}, but couldn't detect {}. exit.".format(col,assembly_dir))
            assembly_col_list.append(col)
            continue

        elif re.match('group\d+$',col):
            group_col+=1
            group_col_list.append(col)
            continue

    if 0 in [samples_col,bam_dir_col,assembly_col,assembly_dir_col,group_col]:
        sys.exit('detected miss components of contig.tsv. exit.')

    if samples_col>1 or bam_dir_col>1:
        sys.exit('detected multiple samples/bam_dir columns. exit.')

    if sum([samples_col,bam_dir_col,assembly_col,assembly_dir_col,group_col]) != len(cols):
        sys.exit('detected unknown column. exit.')

    ## the id of each sample should be unique
    samples=df_config['samples']
    if len(samples)!=len(rm_nan_from_set(set(samples))):
        sys.exit('sample id is not unique. exit.')

    ## the bam of each sample should be unique, and exist
    bams=df_config['bam_dir']
    if len(bams)!=len(rm_nan_from_set(set(bams))):
        sys.exit('bam directory is not unique. exit.')

    ## same assembly should have the same dir, and exist
    for a in assembly_col_list:
        assembly=rm_nan_from_set(set(df_config[a]))
        assembly_file=rm_nan_from_set(set(df_config[a+"_dir"]))

        if len(assembly_file)!=1 or len(assembly)!=1:
            sys.exit("{} has multiple assembly. exit.".format(a))
        if not os.path.exists(next(iter(assembly_file))):
            pass
            #sys.exit("{} doesn't exist. exit.".format(assembly_file))

        assembly_list.append(next(iter(assembly)))
        assembly_loc_list.append(next(iter(assembly_file)))

    ## same group id should retain in same column
    for g in group_col_list:
        group=df_config[g]
        group_set=rm_nan_from_set(set(group))
        for _ in group_set:
            group_name_list.append(_)
            group_list.append(list(df_config.index[df_config[g] == _ ]))
    if len(group_name_list)!=len(rm_nan_from_set(set(group_name_list))):
        sys.exit("group id shows in multiple columns. exit.")

    samples_list=list(df_config['samples'])
    bam_list = list(df_config['bam_dir'])

    return samples_list, assembly_list, assembly_loc_list, group_list, group_name_list,bam_list


