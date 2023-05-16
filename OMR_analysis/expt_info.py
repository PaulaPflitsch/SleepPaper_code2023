#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pyscript.py
#  
#  Copyright 2019 Kumaresh <kumaresh_krishnan@g.harvard.edu>
#  
#  version 1.0

import numpy as np
from hdf5storage import savemat, loadmat
import os, sys

def main():
    
    days = ['2021_06_29']
    fish = [40]
    trials = 20

    control = [0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27,32,33,34,35]
    sleep = [4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,37,38,39]

    #control = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31] # every second fish was put into filtered H2O. odd numbers
    #sleep= [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30] #fish with odd numbers in setup. Here counting starts at 0,
    # so fish with even numbers were stressed with 25mM NaCl

    np.save('expt_info.npy', {'days':days, 'fish':fish, 'trials':trials, 'control':control, 'sleep':sleep})
    
    return 0
    
def createIDMap():
    
    id_map = {}
    
    # Left, right, control
    #id_map['0'] = np.array([-1,50])
    #id_map['1'] = np.array([-1,100])
    #id_map['2'] = np.array([1,50])
    #id_map['3'] = np.array([1,100])
    #id_map['4'] = np.array([0,0])

    #salt: left, right, control
    id_map['0'] = np.array([0, 0])
    id_map['1'] = np.array([-1, 25])
    id_map['2'] = np.array([-1, 50])
    id_map['3'] = np.array([-1, 100])
    id_map['4'] = np.array([0, 0])
    id_map['5'] = np.array([1, 25])
    id_map['6'] = np.array([1, 50])
    id_map['7'] = np.array([1, 100])
    
    # Save ID map
    savemat('ID_map.mat', id_map, format='7.3', oned_as='column', store_python_metadata=True)
    
    return 0


if __name__ == '__main__':

    main()
    createIDMap()
    
    sys.exit()
