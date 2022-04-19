#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 13:48:02 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####

#================================== input =====================================
humann3:
# Gene Family	donkey_3418_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs	donkey_3446_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs	donkey_3611_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs	donkey_4058_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs	donkey_4282_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs	donkey_4687_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs
UNMAPPED	0.921548	0.930192	0.910292	0.914269	0.915457	0.928255
UNGROUPED	0.07377425704386407	0.06668929678790783	0.08420391623811686	0.08154017701510215	0.07942551942900009	0.06775615834239965
K00002	1.88165e-07	3.07085e-07	3.42713e-07	1.66582e-07	6.59669e-08	2.12661e-07

rela_abun_table:
Bacteria,Eukaryota,Archaea
0.968261176998728,0.01592529073003537,0.013180790239979498
0.9783958135574876,0.010658936481966272,0.00865624992670442

#================================== output ====================================

#================================ parameters ==================================
rmUMAPPED
rmUNGROUPED
top=20 # 0 means select all, 20 means select top rich 20

#================================== example ===================================

#================================== warning ===================================
####=======================================================================####
"""
import pandas as pd
import sys

def humann2rel_abun(dp, output,rmUMAPPED = True,rmUNGROUPED = True, top = 0):
    """
    top = 0 means select all, = 20 means select top rich 20
    """
    df_hum = pd.read_csv(dp,sep='\t',index_col=0) # hum is humann

    new_column = [x.replace('_Abundance-RPKs','') for x in df_hum.columns]
    df_hum.columns = new_column
    new_index = [x.replace(',','_') for x in df_hum.index]
    df_hum.index = new_index

    if rmUMAPPED==True:
        df_hum.drop('UNMAPPED',inplace=True,axis=0)
    if rmUNGROUPED==True:
        df_hum.drop('UNGROUPED',inplace=True,axis=0)

    df_hum.T.to_csv(output,sep=',',index = True)

if __name__ == '__main__':
    dp_list = sys.argv[1].split(',')
    output_list = sys.argv[2].split(',')

    for dp,output in zip(dp_list,output_list):
        humann2rel_abun(dp,output)





















