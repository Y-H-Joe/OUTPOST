#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 16:30:24 2022

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
import os
import time
import subprocess
import OUTPOST.check_config
import OUTPOST.decorators
import OUTPOST.functools
import OUTPOST.check_snakefile_config

sys.path.append(os.getcwd())

def wait_until_file_exists(dp_list, error_log, max_circle = 30):
    print(f"OUTPOST: waiting for {dp_list}")
    circle = 0
    for dp in dp_list:
        while circle <= max_circle:
            if os.path.exists(dp):
                time.sleep(5)
                break
            else:
                time.sleep(10)
                circle += 1
            if os.path.exists(error_log):
                return False
            if circle == max_circle:
                print(f"OUTPOST: reached max_circle {max_circle}. quit waiting.")
    return True

def index_genome(genome, bwa, samtools, need_fai = False):
    if not os.path.exists(f"{genome}.bwt.2bit.64"):
        subprocess.run([bwa, 'index', genome])
    if need_fai:
        if not os.path.exists(f"{genome}.fai"):
            subprocess.run([samtools, 'faidx', genome])
            if not os.path.exists(f"{genome}.fai"):
                subprocess.run(['gzip -d', genome])
                subprocess.run(['bgzip', genome])
                subprocess.run([samtools, 'faidx', genome])

def func_begin(*args, **kwargs):
    from OUTPOST.functools import func_begin as fb
    return fb(*args, **kwargs)

def func_end(*args, **kwargs):
    from OUTPOST.functools import func_end as fe
    return fe(*args, **kwargs)
