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

'donkey_3418_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs'
clean_hum_head=2  means only keep donkey_3418

clean_hum_head=2 # 0 means keep all

#================================== example ===================================

#================================== warning ===================================
since I'm not sure what the order of humann3 output table column will be, I seperate
the index file and value file. each of them will be passed to heatmap R plot
####=======================================================================####
"""
import pandas as pd
import re
import sys

def humann2rel_abun(dp):
    #dp=r'../metabolism/humann3/horsedonkey_genefamilies_uniref90names_cpm_{}_unstratified.named.tsv'.format(i)
    dp = r'C:\CurrentProjects\equids_MHC\Prj1\humann3\horsemuledonkeyhinny_genefamilies_uniref90names_cpm_eggnog_unstratified.named.tsv'.format(i)
    #dp_config='config.tsv'
    dp_config = r'C:\CurrentProjects\equids_MHC\Prj1\humann3\equids.prj1.humann3.index'

    rmUMAPPED = True
    rmUNGROUPED = True

    top = 0 # 0 means select all, 20 means select top rich 20
    """
    'donkey_3418_meta_nonEquCab3donkeyviralrRNA_merger1r2_Abundance-RPKs'
    clean_hum_head=2  means only keep donkey_3418
    """
    clean_hum_head=0 # 0 means keep all

    df_hum = pd.read_csv(dp,sep='\t',index_col=0)
    df_config = pd.read_csv(dp_config,sep='\t')

    df_hum_columns = list(df_hum.columns)
    samples=list(df_config['samples'])

    reordered_df_hum_columns=[]
    for s in samples:
        for c in df_hum_columns:
            if bool(re.match(r"{}_?".format(s),c)):
                reordered_df_hum_columns.append(c)
    if len(reordered_df_hum_columns)!=len(samples):
        print("[{}] error. the sample name cannot be reordered. script confuses about the name. exit.".format(sys.argv[0]))
        sys.exit()
    df_hum=df_hum[reordered_df_hum_columns]

    if rmUMAPPED==True:
        df_hum.drop('UNMAPPED',inplace=True,axis=0)
    if rmUNGROUPED==True:
        df_hum.drop('UNGROUPED',inplace=True,axis=0)

    if top!=0:
        suffix='top{}.'.format(str(top))
        df_hum['sum']=df_hum.sum(axis=1)
        df_hum.sort_values(by='sum',inplace=True,ascending=False)
        if top<df_hum.shape[0]:
            df_hum=df_hum.loc[df_hum.index[:top]]
    else:
        suffix=''


    df_hum=df_hum.T
    try:
        df_hum.drop('sum',axis=0,inplace=True)
    except:
        pass

    if clean_hum_head!=0:
        import re
        clean_df_hum_index=['_'.join(re.split('_|-| |\.',x)[:clean_hum_head]) for x in df_hum.index]

    output_value="{}.rel_abun.{}csv".format(dp,suffix)
    df_hum.to_csv(output_value,sep=',',index=True)

if __name__ == '__main__':
    dp_list = sys.argv[1].split(',')
    output_list = sys.argv[2].split(',')

    for dp in zip(dp_list,output_list):
        humann2rel_abun(dp,output)





















