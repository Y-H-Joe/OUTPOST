# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:02:15 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
calculate the relative abundance and counts per million reads for each contig
in each sample.
#================================== input =====================================
#================================== output ====================================
#================================ parameters ==================================
#================================== example ===================================
#================================== warning ===================================
####=======================================================================####
"""
import pandas as pd
import sys

def process_contig_table(dp, samples, output):
    df = pd.read_csv(dp,sep = '\t',index_col = 0)

    sum_ = df[samples].sum()

    def get_rel_abun(numerator, denominator):
        return numerator / denominator

    def get_cpm(numerator, denominator):
        return numerator / denominator * 1000000

    for sample in samples:
        denominator = sum_[sample]
        df[sample + '_rel_abun'] = df[sample].apply(lambda x: get_rel_abun(x, denominator = denominator))
        df[sample + '_cpm'] = df[sample].apply(lambda x: get_cpm(x, denominator = denominator))

    df.to_csv(output, sep = '\t')

if __name__ == '__main__':
    dp = sys.argv[1]
    samples = sys.argv[2].split(',')
    output = sys.argv[3]

    process_contig_table(dp, samples, output)

