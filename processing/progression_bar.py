# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 13:24:57 2018

@author: Simon
"""
from __future__ import print_function
import sys


class my_progression_bar:
    
    def __init__(self, name, total_iteration, h=0.02):
        self.name = name
        self.total_iteration = total_iteration
        self.current_iteration = 0
        self.current_step = 0
        self.h = h
        
    def increment(self):
        assert self.current_iteration <= 1
        
        self.current_iteration += 1./self.total_iteration
        
        if self.current_iteration > self.current_step:
            bar = int(self.current_step/self.h)
            #sys.stdout.write("\r" +" %s |"%self.name+"#"*(bar)+" "*(100-bar)+"| %s pc"%(self.current_step*100))
            print("\r|"+"#"*(bar)+" "*(100-bar)+"|", end='')
            #sys.stdout.flush()
            self.current_step += self.h