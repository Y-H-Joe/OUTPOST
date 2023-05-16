#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 16:12:29 2021

@author: Yihang Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####

#================================== input =====================================

#================================== output ====================================

#================================ parameters ==================================

#================================== example ===================================

#================================== warning ===================================

####=======================================================================####
"""
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import pandas as pd

def res2plot(dp,LDA_cutoff,output):
    df = pd.read_csv(dp,header = None,sep='\t')

    LDA_cutoff = LDA_cutoff

    df1 = df.loc[df[2].notnull()]
    df2 = df1.loc[df1[3] > LDA_cutoff]
    df3 = df2.sort_values(by = [2,3],ascending=False)

    groups=list(set(df3[2]))


    # Create bars
    # barWidth = 0.9
    bars1 = df3.loc[df[2] == groups[0],3].to_list()
    bars2 = df3.loc[df[2] == groups[1],3].to_list()
    # bars1_round=[round(x,2) for x in bars1]
    # bars2_round=[round(x,2) for x in bars2]

    # bars3 = bars1_round + bars2_round

    # The X position of bars
    r1 = list(range(1,len(bars1)+1))
    r2 = list(range(len(bars1)+1,len(bars1)+len(bars2)+1))
    r3 = r1 + r2

    # Create barplot
    left = min(round(min(bars1+bars2),1),0) # set 0 to LDA_cutoff to set the base point

    plt.figure(num=1,figsize=(16,int(df3.shape[0]/5)))
    plt.barh(r1, [x-left for x in bars1], color = 'blue',left =left, label=groups[0])
    plt.barh(r2, [x-left for x in bars2], color = 'yellow',left =left, label=groups[1])
    # Note: the barplot could be created easily. See the barplot section for other examples.

    # Create legend
    plt.legend()
    labels = df3.loc[df[2] == groups[0],0].to_list() + df3.loc[df[2] == groups[1],0].to_list()
    plt.yticks(r3,labels = labels)

    """
    # Text below each barplot with a rotation at 90°
    plt.xticks([r + barWidth for r in range(len(r3))], df3[0].to_list(), rotation=90)

    # Text on the top of each bar
    for i in range(len(r3)):
        plt.text(x = r3[i]-0.2 , y = bars3[i]+0.1, s =bars3[i], size = 6)

    # Adjust the margins
    plt.subplots_adjust()
    """
    # Show graphic
    # plt.show()
    plt.tight_layout()
    plt.ioff()
    plt.savefig(output)
    plt.close()
    # plt.clf()

if __name__=='__main__':
    dp = sys.argv[1]
    output = sys.argv[2]
    LDA_cutoff = int(sys.argv[3])

    try:
        res2plot(dp,LDA_cutoff,output)
    except:
        print("GEMINI: barplots_from_lefse_res_file: ",dp," has no features.skip.")



