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

def headingAngle(raw_data, stimulus, num_bins):

    start = 'bouts_start_stimulus_%03d'%(stimulus)
    end = 'bouts_end_stimulus_%03d'%(stimulus)

    # Compute differences (convention uses start-end)
    all_angles = raw_data[start]['fish_accumulated_orientation'] - \
                       raw_data[end]['fish_accumulated_orientation']

    if all_angles.size == 0:
        return np.array([np.nan]*num_bins)

    # Find bout timestamps and fish location at start
    timestamps = raw_data[start]['timestamp']
    pos_x = raw_data[start]['fish_position_x']
    pos_y = raw_data[start]['fish_position_y']

    # Filter angles based on stimulus time and distance from edge of dish
    vals = np.full(num_bins, np.nan)

    for t in range(0, 20, 2):
        
        lim1, lim2 = t, t+2.0
        locs = np.where((timestamps>lim1) & (timestamps<lim2) & (pos_x**2 + pos_y**2 < 0.81))
        
        angles = all_angles[locs]

        if angles.size == 0:
            continue
        if stimulus == 3:
            vals[int(t//2)] = ((angles < 0).sum() + 0.5*(angles == 0).sum()) / angles.size
        if stimulus == 7:
            vals[int(t//2)] = ((angles > 0).sum() + 0.5*(angles == 0).sum()) / angles.size

    return vals

def extractAngles(experiment,root, num_bins):

    info_path = path.Path() / '..' / experiment

    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days = info['days']
    fish = info['fish']
    trials = info['trials']
    total_fish = np.sum(fish)

    fish_ctr = 0
    stimuli = 2 # Using only 100% left right (3, 7) for this analysis
    num_bins = 10 # Stimulus lasts 20 seconds and bin size of 2 s

    data = np.full((total_fish, trials, stimuli, num_bins), np.nan)

    for day_idx, day in enumerate(days):

        for f in range(fish[day_idx]):
            
            for t in range(trials):

                folder = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp = open(folder, 'rb')
                raw_data = pickle.load(tmp)

                for s_idx, stimulus in enumerate([3,7]):
                    
                    vals = headingAngle(raw_data, stimulus, num_bins)
                    data[fish_ctr, t, s_idx] = vals
                    
                tmp.close()

            fish_ctr += 1
                    
        print(day, fish_ctr, 'fish done')
        
    return data

def processAngles(experiment, data, num_bins):

    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    group_1 = info['control']
    group_2 = info['sleep']

    flat_1 = np.reshape(data[group_1], (-1, num_bins))
    flat_2 = np.reshape(data[group_2], (-1, num_bins))

    avg_data_1 = np.nanmean(flat_1, axis=0)
    avg_data_2 = np.nanmean(flat_2, axis=0)

    sem_data_1 = sem(flat_1, axis=0, nan_policy='omit')
    sem_data_2 = sem(flat_2, axis=0, nan_policy='omit')

    to_save = {}
    to_save['mean_1'] = avg_data_1
    to_save['sem_1'] = sem_data_1.data if np.ma.isMaskedArray(sem_data_1) else sem_data_1
    to_save['mean_2'] = avg_data_2
    to_save['sem_2'] = sem_data_2.data if np.ma.isMaskedArray(sem_data_2) else sem_data_2
    
    return to_save

def plotData(experiment, num_bins):

    data_path = path.Path() / '..' / experiment / f'data_numbouts_timed_{num_bins}'
    tmp = hdf.loadmat(data_path)

    data_1 = tmp['mean_1']
    data_2 = tmp['mean_2']
    sem_1 = tmp['sem_1']
    sem_2 = tmp['sem_2']
    save_dir = path.Path() / '..' / experiment 

    x_range = range(1, num_bins*2, 2)

    os.makedirs(save_dir, exist_ok=True)

    sns.set_style("white")
    sns.set_style('ticks')

    f, ax = plt.subplots()

    ax.errorbar(x_range, data_1, yerr=sem_1, label='control', \
        marker='o', capsize=3.0,color = 'grey')

    ax.errorbar(x_range, data_2, yerr=sem_2, label='light pulses afternoon', \
        marker='o', capsize=3.0, color = 'green')

    ax.axvspan(5,15,0,1, color='grey', alpha=0.1, label='stimulus')

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Avg correctness')
    ax.set_title('Correctness of bouts after stimulus onset (100% coherence) - Time binned')
    ax.set_ylim(0, 1.1)
    ax.set_xlim(0,22)
    ax.legend()
    ax.grid(False)
    sns.despine(top=True, right=True)

    f.savefig(save_dir / f'fig_numbouts_timed_{num_bins}_grey.pdf')
    plt.close(f)

    return 0

def main(experiment, num_bins):

    #root = path.Path() / '..' / '..' / '..' / 'data_hanna_test_06_16_2021' # directory for data
    root = path.Path() / '..' / experiment
    
    data = extractAngles(experiment, root, num_bins)

    to_save = processAngles(experiment, data, num_bins)

    save_dir = path.Path() / '..' / experiment / f'data_numbouts_timed_{num_bins}'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column', store_python_metadata=True)

    return 0

if __name__ == '__main__':

    experiment = 'd8_08_19_2021_recovery'
    num_bins = 10 # num_bins is defined in the extractAngles function anyway for convenience


    main(experiment, num_bins)
    plotData(experiment, num_bins)

    sys.exit(0)
