# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 16:11:19 2012

@author: puya
"""
from numpy import *
import sys

def read_command_line():
    if len(sys.argv) < 5:
        print 'Usage: %s rho_air, v0, T, dt' % \
            sys.argv[0]; sys.exit(1) # abort
    rho_air = float(sys.argv[1])
    v0 = float(sys.argv[2])
    T = float(sys.argv[3])
    dt = float(sys.argv[4])    
    #makeplot = sys.argv[4] in ('on', 'True')
    #dt_values = [float(arg) for arg in sys.argv[5:]]
    return rho_air, v0, T, dt# , makeplot, dt_values



