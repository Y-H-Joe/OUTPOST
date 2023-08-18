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
import time

def print_func_name(func):
    def decorated_func(*args, **kwargs):
        print("{} is running.".format(func.__name__))
        result = func(*args, **kwargs)
        return result
    return decorated_func

def print_run_time(func):  
    def wrapper(*args, **kw):  
        local_time = time.time()  
        value=func(*args, **kw) 
        print('function [%s] run time is %.2f' % (func.__name__ ,time.time() - local_time))
        return value
    return wrapper

