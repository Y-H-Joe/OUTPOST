#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 23:37:40 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####

#================================== input =====================================
streptogramin vat acetyltransferase,glycopeptide resistance gene cluster,AAC(6'),
3.670558122,7.761240366,0.420124122,151.7974789,2.719750897,
#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================

#================================== warning ===================================

####=======================================================================####
"""
import sys
import pandas as pd
import numpy as np
import heapq
import os

def scal_rel_abun(dp,output, rmOthers = False):
    #basename=os.path.basename(dp)
    #output_name=os.path.join(r"../taxa_abun/figs/",basename)
    #dp=r"C:\CurrentProjects\CPS_micro_ts\sample140_abun_top.20.7.rmU.euk.tsv"
    #dp=r"D:\CurrentProjects\equids_MHC\Prj1\humann3\pvalue_matrix_ec_upper_rxn_down.csv"


    ### need to care about the index_col
    df=pd.read_csv(dp,sep=",",index_col=0)

    if rmOthers==True:
        try:
            df.drop('Others',inplace=True,axis=1)
        except:
            rmOthers = False
    ## the global minimum
    second_min=heapq.nsmallest(2,set(df.to_numpy().flatten()))[1]

    ## clean original data
    df.dropna(axis=1, how='all',inplace=True)
    df.replace(0,second_min/100,inplace=True)
    
    ## norm
    # df_norm = (df - df.min()) / (df.max() - df.min())

    ## log10
    df_log10=df.apply(np.log)
    df_log10.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_log10.to_csv(output,index = True,sep = ',')


    ## log2
    # df_log2=df.apply(np.log2)
    # df_log2.replace([np.inf, -np.inf], np.nan, inplace=True)


# =============================================================================
#     if dropna_or_fillmin == 0:
#         suffix1='.dropna'
#     else:
#         suffix1='.fillmin'
#
#     if rmOthers==True:
#         suffix2='.rmOthers'
#     else:
#         suffix2=''
# =============================================================================

    # output_log10=str(dp.strip('.csv')+suffix1+suffix2+".log10"+".csv")

    """
    sns.set(color_codes=True)
    iris = sns.load_dataset("iris")
    species = iris.pop("species")
    """
    """
    ## draw heatmap
    fig,ax = plt.subplots(figsize=(25,5))
    sns.heatmap(df,vmin=0, vmax=1,square=False,cmap="YlGnBu",annot=True)
    fig.savefig(str(dp+".norm.heatmap.pdf"),bbox_inches='tight')

    #设置图片大小
    g= sns.clustermap(df, fmt="d",cmap='YlGnBu')
    ax = g.ax_heatmap
    label_y = ax.get_yticklabels()
    plt.setp(label_y, rotation=360, horizontalalignment='left')
    #设置图片名称，分辨率，并保存
    #plt.savefig('cluster.tif',dpi = 300)
    plt.show()
    """

if __name__=='__main__':
    dp_list = sys.argv[1].split(',')
    # dp_list = [r'cat.rel_abun.normal_vs_obese.at_class.rel_abun.unequal.top20.csv']

    output_list = sys.argv[2].split(',')
    #output_list = [r'cat.rel_abun.normal_vs_obese.at_class.rel_abun.unequal.top20.aa.csv']

    for dp,output in zip(dp_list,output_list):
        #dp=str(r"../taxa_abun/rel_abun/sample12_rel_abun."+str(i)+".rmU.euk.csv.top30.csv")
        #dp=r"../taxa_abun/utest/sample12_rel_abun.{}.rmU.euk.csv_relative_abun_unequal_horse_vs_donkey.csv.top30.csv".format(str(i))
        assert os.path.exists(dp), f"GEMINI: scale_rel_abun_table: {dp} doesn't exist. exit."
        try:
            scal_rel_abun(dp, output)
        except:
            os.system(f"touch {output}")

