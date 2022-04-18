#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 21:22:48 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####

=================================== input =====================================
counts table:
"hmdh_000000000002"	66219	"Bacteria"	"Firmicutes"	"Clostridia"	"Eubacteriales"	"Lachnospiraceae"	"Lachnoclostridium"	"Lachnoclostridium phytofermentans"	0	0	1	0	0	0	0	0	0	2	0	0	0	0	0	0	2	3	8
"hmdh_000000000003"	1265	"Bacteria"	"Firmicutes"	"Clostridia"	"Eubacteriales"	"Oscillospiraceae"	"Ruminococcus"	"Ruminococcus flavefaciens"	2	0	0	0	0	0	0	0	0	1	1	0	0	0	3	0	1	0	8
"hmdh_000000000004"	0	""	""	""	""	""	""	""	0	0	0	1	0	0	0	0	0	0	0	0	0	0	2	0	1	0	4
=================================== output ====================================

================================= parameters ==================================

=================================== example ===================================

=================================== warning ===================================

####=======================================================================####
"""
import pandas as pd
import sys
import os
from scipy import stats

def pvalue(pair,list1,list2):
    if pair==True:
        k,p=stats.wilcoxon(list1,list2,alternative="two-sided")
    else:
        k,p=stats.mannwhitneyu(list1,list2,alternative="two-sided")
    return p

def p_adjust_bh(p):
    import numpy as np
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):
    """
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])
    """
    from numpy import array, empty
    pvalues = array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in range(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
    return new_pvalues


def list_slice(list1,list2):
    return [list1[x] for x in list2]

if __name__ == '__main__':
    try:
        ## parameters
        dp = sys.argv[1]
        group1 = [int(x) for x in sys.argv[2].split(',')]
        group2 = [int(x) for x in sys.argv[3].split(',')]
        output = sys.argv[4]
        sample_list = sys.argv[5].split(',')

        sample_size=len(sample_list)

        count_col=range(9,9+sample_size)

        pair=False
        filter_ = False

        ## begin
        df=pd.read_csv(dp,sep='\t',header=None)

        group1_col=list_slice(count_col, group1)
        group2_col=list_slice(count_col, group2)


        p_list=[]
        for i in df.index:
            list1=list(df.loc[i,group1_col])
            list2=list(df.loc[i,group2_col])
            p_list.append(pvalue(pair, list1, list2))

        #adj_p_list=p_adjust_bh(p_list)
        adj_p_list=correct_pvalues_for_multiple_testing(p_list)
        df.columns = ['contigID','taxaID','superkingdom','phylum','class','order','family','genus','species'] + \
                        sample_list
        df['pvalue']=p_list
        df['qvalue']=adj_p_list

        ###############################################################################
        ## filter

        # remove taxa without family annotation
        # remove pvalue > 0.1
        # remove 0 counts number > sample_size-3
        # remove qvalue > 0.2
        if filter_:
            left_index=[]
            for i in df.index:
                ## family col is 6
                ## sum col is 27
                if df.loc[i].isnull()[6] or df.loc[i,'pvalue']>0.1 or sum(df.loc[i,count_col]==0)>(sample_size-3) or df.loc[i,'qvalue'] > 0.2:
                    continue
                else:
                    left_index.append(i)
            df_left=df.loc[left_index]
            df_left.to_csv(output,sep='\t',index=None)
        else: df.to_csv(output,sep='\t',index=None)
    except Exception as e:
        import traceback
        error_log = sys.argv[6]
        os.system("touch " + error_log)
        print(f"GEMINI: {e}")
        traceback.print_exc()



