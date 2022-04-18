# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 22:14:46 2020

@author: Y.H. Zhou

sample16_head600.tsv

Bosmeta24_000000000010	Bacteria	Firmicutes	Clostridia	NA	NA	NA	Clostridia bacterium	4	0	0	0	0	3	1	2	5	3	1	0	0	2	1	1	23
Bosmeta24_000000000015								2	0	1	0	0	0	1	0	3	2	2	0	3	2	1	0	17
Bosmeta24_000000000017	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Sarcina	Sarcina sp. DSM 11001	0	6	0	0	1	2	0	0	0	0	0	2	0	1	0	5	17
Bosmeta24_000000000018	Bacteria	Firmicutes	Bacilli	Bacillales	Paenibacillaceae	Paenibacillus	Paenibacillus albus	2	0	0	1	5	2	0	0	0	0	0	2	0	0	0	2	14
"""

import pandas as pd
import numpy as np
import time
import math
import sys

def print_run_time(func):  
    def wrapper(*args, **kw):  
        local_time = time.time()  
        value=func(*args, **kw) 
        print('function [%s] run time is %.2f' % (func.__name__ ,time.time() - local_time))
        return value
    return wrapper

#@print_run_time
def df_rman(dp):
    df=pd.read_csv(dp,header=None,index_col=None,sep="\t")
    df_rmnan=df.replace(np.nan,"Unknown" )
    return df_rmnan

#@print_run_time
def taxa_count_table_to_sorted_list_dict(df_rmnan,sample_N,taxo_level_N,level):
    #k_set=[]
    #for k in range(1,1+taxo_level_N): # taxo levels
    #    k_set.append(set(list(df_rmnan[k]))) ## redundanced taxas
      
    ## relative abundance counts in terms of sample and taxo
    df_rmnan_groupby=df_rmnan.loc[:,[level]+list(range(taxo_level_N+1,taxo_level_N+1+sample_N))].groupby([level]).sum()
    taxa_count_pairs=df_rmnan_groupby.to_dict()
    samples=[taxa_count_pairs[c] for c in taxa_count_pairs.keys() ]
    #k_counts=samples[-1]
    #del(samples[-1])
        
    # sort sample counts and get abundance
    samples_sum=[sum(s.values()) for s in samples]
    samples_abun_sort=[]
    for s in range(sample_N):
        sample_abun={}
        sample=samples[s]
        for key in sample.keys():
            if type(key) is not str:
                #print("occured a nan key:",key)
                if math.isnan(key):
                    print("in ",s+1,"th"," sample.")
                    print("sample[key]: ",sample[key])
                    print("samples_sum[s]: ",samples_sum[s])
            else:
                pass
            try:
                sample_abun[key]=float(int(sample[key])/int(samples_sum[s]))
            except:
                print("cannot calculate abundance.")
                print("in ",s,"th"," sample.")
                print("key: ",key)
                print("sample[key]: ",sample[key])
                print("samples_sum[s]: ",samples_sum[s])
        sample_abun_sor={k: v for k, v in sorted(sample_abun.items(), key=lambda item: item[1],reverse=True)}
        samples_abun_sort.append(sample_abun_sor) # sorted, large to samll
        
    return samples_abun_sort


#@print_run_time
def rmU_adj_samples_abun_sort(samples_abun_sort,taxo):
    samples_abun_rmU=[]
    for d in samples_abun_sort:
        d_Copy=d.copy()
        try:
            d_Copy.pop(taxo)
        except:
            pass
        samples_abun_rmU.append(d_Copy)
    
    samples_abun_sort_rmU_adj=[] #adjusted
    for sample in samples_abun_rmU:
        sample_adj={}
        abun_sum=sum(sample.values()) #0.78156....
        for key in sample.keys():
            sample_adj[key]=float(sample[key]/abun_sum)
        samples_abun_sort_rmU_adj.append(sample_adj)
    return samples_abun_sort_rmU_adj


#@print_run_time
def extract_samples_abun_top(samples_abun_sort,top):
    samples_abun_top=[]
    for s in range(len(samples_abun_sort)):
        s_abun={}
        sample=samples_abun_sort[s]
        top_taxo=list(sample.keys())[:top-1]
        for i in range(top-1):
            taxo=top_taxo[i]
            s_abun[taxo]=sample[taxo]
        samples_abun_top.append(s_abun)
    
    return(samples_abun_top)


# top mode
#@print_run_time
def df_samples_abun_top(samples_abun_sort,top=0):
    # top mode / all mode
    """
    start=time.time()
    if top>0:
        samples_abun_top=extract_samples_abun_top(samples_abun_sort,top)
    elif top==0:
        All=len(samples_abun_sort[0])
        samples_abun_top=extract_samples_abun_top(samples_abun_sort,All)
    
    print(time.time()-start)
    
    start=time.time()
    samples_abun_top_df=pd.DataFrame.from_dict(samples_abun_top)
    index=list(samples_abun_top_df.index)
    col=list(samples_abun_top_df.columns)
    for i in index: # make up nan cause samples have various top10
        for c in col:
            if np.isnan(samples_abun_top_df.loc[i,c]):
                samples_abun_top_df.loc[i,c]=samples_abun_sort[i][c]
    print(time.time()-start)
    
    start=time.time()
    if top>0: # top mode
        samples_abun_top_df["Others"]=0
        for i in index:
            samples_abun_top_df.loc[i,"Others"]=float(1-sum(samples_abun_top_df.loc[i]))
    print(time.time()-start)
    samples_abun_top_df.columns=["Unknown" if taxo !=taxo else taxo for taxo in samples_abun_top_df.columns] # nan != nan trick
    """
    samples_abun_top=samples_abun_sort
    samples_abun_top_df=pd.DataFrame.from_dict(samples_abun_top)
    samples_abun_top_df.columns=["Unknown" if taxo !=taxo else taxo for taxo in samples_abun_top_df.columns]
    
    return samples_abun_top_df

if __name__=='__main__':
    dp=r"hd.align2contigs.hmdh.counts.tsv"
    dp_config="config.tsv"
    
    rmUnknown=True
    kaiju="euk" # "nr","euk"
    taxo_level_N=8 ## total levels of taxonomy
    sample_N=12
    
    df_rmnan=df_rman(dp)
    df_config=pd.read_csv(dp_config,sep='\t',index_col=0)
    for level in range(1,taxo_level_N+1):
        samples_abun_sort=taxa_count_table_to_sorted_list_dict(df_rmnan,sample_N,taxo_level_N,level)
        if rmUnknown!=True:
            samples_abun_top_df=df_samples_abun_top(samples_abun_sort)
        if rmUnknown==True:
            samples_abun_sort_rmU_adj=rmU_adj_samples_abun_sort(samples_abun_sort,"Unknown")   
            samples_abun_top_df=df_samples_abun_top(samples_abun_sort_rmU_adj)
        try:
            samples_abun_top_df.index=df_config.index
            samples_abun_top_df.to_csv(str("sample"+str(sample_N)+"_rel_abun."+str(level)+[".",".rmU."][bool(rmUnknown)]+kaiju+".csv"),sep=",",header=True,index=True) 
        except:
            import sys
            print("config error.exit.")
            sys.exit()
    


