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
import re
import os
import time
import OUTPOST.check_config
import OUTPOST.decorators
import OUTPOST.functools

sys.path.append(os.getcwd())
check_config = OUTPOST.check_config.check_config
func_begin = OUTPOST.functools.func_begin
func_end = OUTPOST.functools.func_end

def wait_unti_file_exists(dp_list, error_log, max_circle = 30):
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
