#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 14:40:41 2021

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
import pandas as pd
import os

if not os.path.exists('../taxa_abun/alpha_beta_diversity/'):
    os.makedirs('../taxa_abun/alpha_beta_diversity/')

dp_config="config.tsv"
dp_rel_abun="../taxa_abun/rel_abun/sample12_rel_abun.8.rmU.euk.csv"
rel_table_basename=os.path.basename(dp_rel_abun)

output_group="alpha_beta_diversity_group.tsv"
output_rel_abun="../taxa_abun/alpha_beta_diversity/"+rel_table_basename+".alpha_beta.csv"

df_config=pd.read_csv(dp_config,sep='\t')
df_rel_abun=pd.read_csv(dp_rel_abun)

group_col=[x for x in df_config.columns if 'group' in x ]

df_new_rel_table=pd.DataFrame()
df_group=pd.DataFrame()
for g in group_col:
    index_of_g=df_config.index[df_config[g].notna()]
    df_new_rel_table=pd.concat([df_new_rel_table,df_rel_abun.loc[index_of_g]])
    
    df_contig_tmp=df_config.loc[index_of_g]
    df_group_tmp=pd.DataFrame()
    df_group_tmp['samples']=df_contig_tmp[['samples']]
    df_group_tmp['group']=df_contig_tmp[[g]]
    df_group_tmp['others']=['A']*df_group_tmp.shape[0]
    df_group=pd.concat([df_group,df_group_tmp])
    

df_group.index=range(df_group.shape[0])
# to make sure each sample is unique in alpha-beta diversity analysis
for i in df_group.index:
    df_group.loc[i,'samples']=str(df_group.loc[i,'samples'])+"_"+str(i)

df_new_rel_table.to_csv(output_rel_abun,index=None)
df_group.to_csv(output_group,sep='\t',index=None)


