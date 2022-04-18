# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 18:21:11 2020

@author: Y.H. Zhou

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ discription ==============================####
#
#================================== input =====================================
U	horsedonkeymeta18_000000000069	0
C	horsedonkeymeta18_000000000191	1972457	106	1972457,	OQY42018.1,	NDEPSKASRPFNADRDGFVMGEG,	Bacteria; Fusobacteria; Fusobacteriia; NA; NA; NA; Fusobacteriia bacterium 4572_74
#================================== output ====================================
hmdh_000000000003	1265	Bacteria	Firmicutes	Clostridia	Eubacteriales	Oscillospiraceae	Ruminococcus	Ruminococcus flavefaciens
hmdh_000000000004	0\t\t\t\t...

#================================ parameters ==================================
#
#================================== example ===================================
#python3 thiscript.py CPS_meta22_contigs.fa.euk.nm
#================================== warning ===================================
#
####=======================================================================####
"""
# import pandas as pd
import sys

infile=sys.argv[1]

with open(infile,'r') as r, open(str(infile+".tsv"),'w') as w:
    for line in r.readlines():
        line0 = line.strip()
        line_list = line0.split('\t')
        contigID = line_list[1]
        taxaID = line_list[2]
        if taxaID == '0':
            taxa = '\t'*6 + '\n'
        else:
            taxa = '\t'.join([x.strip(' ') for x in line_list[-1].strip(';').split(';')]) + '\n'
        w.write(contigID + '\t' + taxaID + '\t' + taxa)



"""
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
"""
