#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:37:53 2022

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
import pandas as pd
import os
import re
import sys
import mwgx_pipeline as mp
config="mwgx_pipeline/config.tsv"
df_config=pd.read_csv(config,sep='\t')

assembly_list, assembly_loc_list, group_list, group_name_list=mp.check_config(df_config)
