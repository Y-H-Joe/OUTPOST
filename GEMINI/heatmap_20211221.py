#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 20:01:41 2021

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
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42


if not os.path.exists('../taxa_abun/figs/'):
    os.makedirs('../taxa_abun/figs/')
    
for i in range(7,8):
    try:
        #dp=str("../taxa_abun/utest/sample12_rel_abun."+str(i)+".rmU.euk.csv_relative_abun_unequal_horse_vs_donkey.csv.top30.csv.log10.fillmin.csv")
        #dp=str("../taxa_abun/rel_abun/sample12_rel_abun."+str(i)+".rmU.euk.csv.top30.csv.log10.fillmin.csv")
        dp = r"D:\CurrentProjects\equids_MHC\Prj1\viruses\sample18_rel_abun.{}.rmU.euk.csv.log10.fillmin.csv.top100.csv".format(i)
        basename=os.path.basename(dp)
        #output=str("../taxa_abun/figs/"+basename+".heatmap.pdf")
        output=r"D:\CurrentProjects\equids_MHC\Prj1\viruses\sample18_rel_abun.{}.rmU.euk.log10.fillmin.top100.heatmap.pdf".format(i)
        
        df=pd.read_csv(dp,index_col=0)
        
        sns.set_theme(color_codes=True)
        
        g = sns.clustermap(df,col_cluster=False,xticklabels=False)
        
        #g.ax_heatmap.tick_params( bottom=False)
        #plt.tight_layout()
        plt.savefig(output)
    except:
        print("Cannot find ",dp,". Skip.")

