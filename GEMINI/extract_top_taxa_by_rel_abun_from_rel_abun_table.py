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
This script is to strictly select 20 most abundant
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
import os
import sys

if  __name__ == '__main__':
    try:
        dp_list = sys.argv[1].split(',')
        output_list = sys.argv[2].split(',')
        top = int(sys.argv[3])

        for dp,output_name in zip(dp_list,output_list):
            print(f"dp: {dp}")
            print(f"output_name: {output_name}")
            try:
                df = pd.read_csv(dp,sep=",",header=0,index_col=0)
            except:
                os.system(f'cp {dp} {output_name}')
                continue
            sum1 = df.sum(axis=0)
            sum1_sort = sum1.sort_values(ascending=False) ## sort here
            # print(f"sum1_sort: {sum1_sort}")
            # print(f"sum1: {sum1}")
            if len(sum1_sort) > top :
                sum1_sort_top = sum1_sort[:top]
                df_top = df[sum1_sort_top.keys()]
                df_top.to_csv(output_name,index=True)
            else:
                print("GEMINI: ",dp," has no enough taxa to extract. remain original file.")
                df.to_csv(output_name,index=True)


    except Exception as e:
        import traceback
        error_log = sys.argv[4]
        os.system("touch " + error_log)
        print(f"GEMINI: {e}")
        traceback.print_exc()

