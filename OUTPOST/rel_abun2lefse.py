# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 17:59:38 2020

@author: Y.H. Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
#
#================================== input =====================================
#
#================================== output ====================================
#
#================================ parameters ==================================
#
#================================== example ===================================
#
#================================== warning ===================================
# seems the length of head can only be 3
####=======================================================================####
"""
import pandas as pd
import sys
import os

def rel_abun2lefse(dp, output, class_):
    df = pd.read_csv(dp,sep=",",index_col=0).T
    subClass = ['subClass']*len(class_)
    Subject = [f"{sample}_{group}" for sample,group in zip(list(df.columns),class_)]
    df.columns=range(df.shape[1])

    head_index = ['Class','subClass','Subject']
    head_data = [class_,subClass,Subject]
    head_dict = dict(zip(head_index,head_data))
    head_df = pd.DataFrame.from_dict(head_dict,orient='index')

    lefse_df = pd.concat([head_df,df],axis=0)

    lefse_df.to_csv(output,header=None,sep='\t')


if __name__=='__main__':
    dp_list = sys.argv[1].split(',')
    output_list = sys.argv[2].split(',')
    class_ = sys.argv[3].split(',')

    for dp,output in zip(dp_list,output_list):
        try:
            rel_abun2lefse(dp, output, class_)
        except:
            os.system(f"touch {output}")


