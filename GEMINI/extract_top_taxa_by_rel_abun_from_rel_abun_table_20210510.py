#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 16:27:54 2021

@author: Y.H. Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
#
#================================== input =====================================
rel abun table-
,Spirotrichea,Rhodothermia
donkey1,0.0,0.0
donkey2,0.2,0.0
donkey3,0.0,0.3

#================================== output ====================================
my top XX relative abundance tables are generated selecting union set, therefore the number of
taxa in top 20 is actually larger than 20. This script is to strictly select 20 most abundant
taxa from relative abundance table using sum of relative abundance.
#================================ parameters ==================================
#
#================================== example ===================================
this script automatically sort from largest abun to lowest


#================================== warning ===================================
#
####=======================================================================####
"""
import pandas as pd

for i in range(1,9):
    dp=r"..\taxa_abun\utest\sample12_rel_abun.{}.rmU.euk.csv_relative_abun_unequal_horse_vs_donkey.csv".format(str(i))
    #dp=str(r"../taxa_abun/rel_abun/sample12_rel_abun."+str(i)+".rmU.euk.csv")

    top=30

    df=pd.read_csv(dp,sep=",",header=0,index_col=0)

    sum1=df.sum(axis=0)
    sum1_sort=sum1.sort_values(ascending=False) ## sort here

    if len(sum1_sort) > top :
        sum1_sort_top=sum1_sort[:top]
        output_name=str(dp+".top"+str(top)+".csv")
        df_top=df[sum1_sort_top.keys()]
        df_top.to_csv(output_name,index=True)
    else:
        print(dp," has no enough taxa to extract. Skip.")


