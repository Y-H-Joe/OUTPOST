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
import GEMINI.check_config
import GEMINI.decorators
import GEMINI.functools

sys.path.append(os.getcwd())
check_config = GEMINI.check_config.check_config
func_begin = GEMINI.functools.func_begin
func_end = GEMINI.functools.func_end

def wait_unti_file_exists(dp_list, error_log, max_circle = 1800):
    print(f"GEMINI: waiting for {dp_list}")
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
    return True
    