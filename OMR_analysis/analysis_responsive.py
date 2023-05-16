#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pyscript.py
#  
#  Copyright 2020 Kumaresh <kumaresh_krishnan@g.harvard.edu>
#
#  version 1.0
#  

import os, sys
import numpy as np

import hdf5storage as hdf
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

import path
import pickle

from matplotlib import cm
from colorspacious import cspace_converter

def findResponse(raw_data, stimulus):

    start = 'bouts_start_stimulus_%03d'%(stimulus)
    end = 'bouts_end_stimulus_%03d'%(stimulus)

    # Did they respond at all - then there should be bouts
    angles = raw_data[start]['fish_accumulated_orientation']
    
    if angles.size == 0:
        return 0

    return 1

def extractAngles(experiment,root):

    info_path = path.Path() / '..' / experiment

    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days = info['days']
    fish = info['fish']
    trials = info['trials']
    total_fish = np.sum(fish)

    fish_ctr = 0
    stimuli = 8

    data = np.zeros((total_fish, trials, stimuli))

    for day_idx, day in enumerate(days):

        for f in range(fish[day_idx]):
            
            for t in range(trials):

                folder = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp = open(folder, 'rb')
                raw_data = pickle.load(tmp)

                for stimulus in range(stimuli):
                    
                    response = findResponse(raw_data, stimulus)
                    data[fish_ctr, t, stimulus] = response
                    
                tmp.close()

            fish_ctr += 1
                    
        print(day, fish_ctr, 'fish done')
        
    return data

def processData(experiment, data):

    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    group_1 = info['control']
    group_2 = info['sleep']

    data_1 = data[group_1]
    data_2 = data[group_2]

    responses_1 = data_1.sum(axis=1) / data_1.shape[1]
    responses_2 = data_2.sum(axis=1) / data_2.shape[1]

    to_save = {}
    to_save['responses_1'] = responses_1
    to_save['responses_2'] = responses_2

    return to_save

def main(experiment):

    #root = path.Path() / '..' / '..' / '..' / 'data_hanna_test_06_16_2021' # directory for data
    root = path.Path() / '..' / experiment
    
    data = extractAngles(experiment, root)

    to_save = processData(experiment, data)

    save_dir = path.Path() / '..' / experiment / f'data_responses'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column', store_python_metadata=True)

    return 0

def plotResults(experiment):

    data_path = path.Path() / '..' / experiment / 'data_responses.mat'
    save_dir = path.Path() / '..' / experiment

    tmp = hdf.loadmat(data_path)
    data_1 = tmp['responses_1']
    data_2 = tmp['responses_2']

    sns.set_style('white')
    sns.set_style('ticks')
    
    f, (ax1, ax2) = plt.subplots(1,2, figsize=(10,12))

    m1 = ax1.imshow(data_1, cmap='Reds', vmin=0, vmax=1.0, label='control')
    m2 = ax2.imshow(data_2, cmap='Reds', vmin=0, vmax=1.0, label='sleep deprived')
    plt.suptitle('Percentage response of fish for each stimulus - 30 trials')
    f.colorbar(m1, ax=ax1); f.colorbar(m2, ax=ax2)
    ax1.set_title(f'Control fish: {data_1.mean() * 100:.2f}%')
    ax2.set_title(f'Sleep deprived fish (light pulses afternoon): {data_2.mean() * 100:.2f}%')
    ax1.set_xlabel('Stimulus'); ax2.set_xlabel('Stimulus')
    ax1.set_ylabel('Fish number'); ax2.set_ylabel('Fish number')

    f.savefig(save_dir / 'response_percentage.pdf')
    plt.close(f)

    return 0

if __name__ == '__main__':

    experiment = 'd8_08_19_2021_recovery'

    main(experiment)
    plotResults(experiment)

    sys.exit(0)
