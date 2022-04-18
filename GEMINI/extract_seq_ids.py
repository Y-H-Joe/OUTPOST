# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:30:35 2020

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
# awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' head100.fq > head100.ids.len
####=======================================================================####
"""
import sys
dp=sys.argv[1]
out=sys.argv[2]

contig_info=True
"""
>Bosmeta24_000000000010
479
>B...
"""


start="1"
end=0
with open(dp,'r') as f:
    with open(out,'w') as new_f:
        lines=f.readlines()
        for line in lines:
            if line[0]==">":
                new_line=str(line[1:]).strip()
                new_f.write(new_line)
                new_f.write('\t')
                new_f.write(start)
                new_f.write('\t')
            else:
                end=str(line).strip()
                new_f.write(end)
                new_f.write('\n')


                    
            