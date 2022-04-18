# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:21:11 2020

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
#python3 thiscript.py CPS_meta22_contigs.fa.euk.nm 
#================================== warning ===================================
#
####=======================================================================####
"""
import pandas as pd
import sys

infile=sys.argv[1]
df=pd.read_csv(infile,header=None,index_col=None)
index=df.index
rmtaxaID=False
lis_tap=[]
for i in index:
    i1=list(df.loc[i])
    i2=i1[0].strip().split(';')
    i3=[x.strip().split('\t') for x in i2]
    i4=[x for y in i3 for x in y]
    del(i4[0]) ## remove the U/C column
    if rmtaxaID:
        del(i4[1]) ## remove the taxon ID column
    lis_tap.append(i4)

df_tap=pd.DataFrame(lis_tap)
if rmtaxaID:
    df_tap.drop([8],axis=1,inplace=True)
else:
    df_tap.drop([9],axis=1,inplace=True)
df_tap.to_csv(str(infile+".tsv"),header=None,index=None,sep='\t')

