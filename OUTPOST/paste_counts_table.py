#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 15:17:25 2020

@author: 23712

Contact: yihangjoe@foxmail.com
         https://github.com/Y-H-Joe/

####============================ description ==============================####
idx file:
hmdh_nucl_000000001	96	0	0
hmdh_nucl_000000002	156	0	0

the 4 columns are: sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
we need the 3rd column.
#================================== input =====================================
#================================== output ====================================
#================================ parameters ==================================
#================================== example ===================================
#================================== warning ===================================
####=======================================================================####
"""
import os
import sys
import pandas as pd

col = sys.argv[1] ## col usually is 4, because the 4th col is the counts col
names = sys.argv[2].split(',')
output = sys.argv[3]
config = sys.argv[4]

rm_batch_effect = eval(sys.argv[5])

df_config=pd.read_csv(config,sep='\t')
res = {}
output_tmp = output+'_tmp'
output_tmp2 = output + '_has_batch_effect'
## the first one need specially processed to keep the first column CPSmeta22_000000000001
## for further merge with kaiju counts table
os.system(str('cut -f '+"1,"+col+" "+names[0]+" > "+names[0]+".pure"))

for name in names[1:]:
    ## intermediate file name=name+".pure"
    os.system(str('cut -f '+col+" "+name+" > "+name+".pure"))

input_paste=" ".join([str(name+".pure") for name in names])
os.system(str("paste"+" "+input_paste+" > "+output_tmp))


# check output
name_num = len(names)
with open(output_tmp,'r') as r,open(output_tmp2,'w') as w:
    for line in r.readlines():
        items = line.strip().split('\t')
        if len(items)-1 != name_num:
            sys.exit("OUTPOST: the bam files were not mapped to the same assembly, or truncated.")
        if items[0] != '*':
            w.write(line)

os.system(f"rm {output_tmp}")

# update in v1.1
# use CPM to remove the batch effect
# apply CPM to each group
if rm_batch_effect:
    batch_list = list(df_config['batch'])
    batch_num = len(set(batch_list))
    batch_dict = {}
    for idx,batch in enumerate(batch_list):
        if batch not in batch_dict.keys():
            # +1 because the 0 column is 'mammals_000000000002'
            batch_dict[batch] = [idx + 1]
        else:
            batch_dict[batch].append(idx + 1)
    
    df = pd.read_csv(output_tmp2, header = None, sep = '\t', index_col = 0)
    for batch,idxes in batch_dict.items():
        batch_sum = sum(df.loc[:,idxes].sum(axis = 0))
        df.loc[:,idxes] = df.loc[:,idxes].mul(1000000/batch_sum) # CPM
    df.to_csv(output,header = None,sep = '\t', index=True)
else:
    os.system(f"mv {output_tmp2} {output}")

"""
# invoke with column nr to extract as first parameter followed by
# filenames. The files should all have the same number of rows

import sys
import pandas as pd
from contextlib import ExitStack
import fileinput

col = int(sys.argv[1])
namefile=sys.argv[2]
output=sys.argv[3]
res = {}

names_df=pd.read_csv(namefile,header=None)
names=list(names_df[0])
names_tup=tuple(set(names))



with open(output,'w') as o:
    with ExitStack() as stack:
        files = [stack.enter_context(open(name)) for name in names]
        for line_nr,line in enumerate(files[0]):
            output_line='\t'.join([file])



    pass

for file_name in names:
    for line_nr, line in enumerate(open(file_name)):
        res.setdefault(line_nr, []).append(line.strip().split('\t')[col-1])


    for line_nr in sorted(res):
        o.writelines('\t'.join(res[line_nr]))

"""