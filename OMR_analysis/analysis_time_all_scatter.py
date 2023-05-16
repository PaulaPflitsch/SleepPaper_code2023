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

def headingAngle(raw_data, stimulus):

    start = 'bouts_start_stimulus_%03d'%(stimulus)
    end = 'bouts_end_stimulus_%03d'%(stimulus)
    
    # Find bout timestamps and fish location at start
    timestamps = raw_data[start]['timestamp']
    pos_x = raw_data[start]['fish_position_x']
    pos_y = raw_data[start]['fish_position_y']

    angles = raw_data[start]['fish_accumulated_orientation'] - \
                       raw_data[end]['fish_accumulated_orientation']

    if angles.size == 0:
        return np.array([np.nan]*50), np.nan

    # Find time to first bout

    lim1, lim2 = 5.0, 15.00

    filt = (timestamps>lim1) & (timestamps<lim2) & (pos_x**2 + pos_y**2 < 0.81)
    
    bouts = (filt == 1)

    if bouts.sum() == 0:
        return np.array([np.nan]*50), np.nan
    else:
        first = timestamps[np.where(bouts)[0][0]] - lim1
        vals, _ = np.histogram(first, range=(0,3), bins=50)
        vals[0] = 0
        
        return vals, first

def extractAngles(experiment,root):

    info_path = path.Path() / '..' / experiment

    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days = info['days']
    fish = info['fish']
    trials = info['trials']
    total_fish = np.sum(fish)

    fish_ctr = 0
    stimuli = 2 #info['stimuli'] # Only doing 100%

    data_first = np.zeros((total_fish, trials, stimuli, 50))
    data_first_raw = np.zeros((total_fish, trials, stimuli))

    for day_idx, day in enumerate(days):

        for f in range(fish[day_idx]):

            for t in range(trials):

                folder = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp = open(folder, 'rb')
                raw_data = pickle.load(tmp)

                for stimulus in [3,7]: # Only doing 100%
                    
                    first,val = headingAngle(raw_data, stimulus)
                    data_first[fish_ctr, t, stimulus%3] = first
                    data_first_raw[fish_ctr, t, stimulus % 3] = val

                tmp.close()

            fish_ctr += 1
                    
        print(day, fish_ctr, 'fish done')
        
    return data_first, data_first_raw

def processAngles(experiment, data_first, data_first_raw):

    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    group_1 = info['control']
    group_2 = info['sleep']

    doubled = np.nansum(data_first, axis=2)
    
    first_1 = doubled[group_1]
    first_2 = doubled[group_2]

    avg_first_1 = np.nansum(first_1, axis=1)
    avg_first_2 = np.nansum(first_2, axis=1)

    mean_first_1 = np.nanmean(avg_first_1, axis=0)
    mean_first_2 = np.nanmean(avg_first_2, axis=0)

    sem_first_1 = sem(avg_first_1, axis=0, nan_policy='omit')
    sem_first_2 = sem(avg_first_2, axis=0, nan_policy='omit')

    to_save = {}
    to_save['avg_first_1'] = mean_first_1
    to_save['sem_first_1'] = sem_first_1.data if np.ma.isMaskedArray(sem_first_1) else sem_first_1
    to_save['avg_first_2'] = mean_first_2
    to_save['sem_first_2'] = sem_first_2.data if np.ma.isMaskedArray(sem_first_2) else sem_first_2
    to_save['raw_first_1'] = avg_first_1
    to_save['raw_first_2'] = avg_first_2
    to_save['raw_value_1'] = data_first_raw[group_1]   #average across stimuli and trials and take the mean
    to_save['raw_value_2'] = data_first_raw[group_2]    # average across stimuli and trials and take the mean

    
    return to_save


# https://stackoverflow.com/questions/60229586/how-to-make-jitterplot-on-matplolib-python
def jitter_dots(dots):
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.005, 0.005, offsets.shape[0])
    dots.set_offsets(jittered_offsets)

def main(experiment):

    #root = path.Path() / '..' / '..' / '..' / 'data_hanna_test_06_16_2021' # directory for data
    #root = path.Path() / '..' / '..' / '..' / 'data_combined_d8'
    root = path.Path() / '..' / experiment

    data_first, data_first_raw = extractAngles(experiment, root)

    to_save = processAngles(experiment, data_first, data_first_raw)

    save_dir = path.Path() / '..' / experiment / f'data_time'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column', store_python_metadata=True)

    return 0

def plotHistogram(experiment, prob=False):

    data_path = path.Path() / '..' / experiment / f'data_time'
    tmp = hdf.loadmat(data_path)

    first_1 = tmp['avg_first_1']
    first_2 = tmp['avg_first_2']

    norm_1 = first_1.sum()
    norm_2 = first_2.sum()

    first_1 /= norm_1
    first_2 /= norm_2

    sem_first_1 = tmp['sem_first_1'] / norm_1
    sem_first_2 = tmp['sem_first_2'] / norm_2

    load_raw_1 = tmp['raw_value_1']
    load_raw_2 = tmp['raw_value_2']

    raw_1 = np.nanmean(load_raw_1, axis=(1,2))
    raw_2 = np.nanmean(load_raw_2, axis=(1,2))

    save_dir = path.Path() / '..' / experiment

    x_vals = np.linspace(0,3000,50)

    x_range = range(first_1.shape[0])

    sns.set_style('white')
    sns.set_style('ticks')

    f, ax = plt.subplots()

    ax.plot(x_vals, first_1, label='control', color='black')
    ax.plot(x_vals, first_2, label='0,1µM Melatonin', color='orange')


    plt.fill_between(x_vals, first_1-sem_first_1, first_1+sem_first_1, alpha=0.3, color='grey')
    plt.fill_between(x_vals, first_2-sem_first_2, first_2+sem_first_2, alpha=0.3, color='orange')

    ax.set_xlabel(f'Time (ms)')
    ax.set_ylabel(f'Frequency')
    ax.set_title(f'Time to first bout')
    ax.set_xlim(0,3000)
    ax.set_ylim(0,0.10)
    
    ax.legend()
    ax.grid(False)
    sns.despine(ax=ax, top=True, right=True)

    f.savefig(save_dir / f'fig_first.pdf')

    plt.close(f)

    mean_1 = np.nanmean(raw_1*1000)
    mean_sem_1 = sem(raw_1*1000, nan_policy='omit')
    mean_2 = np.nanmean(raw_2*1000)
    mean_sem_2 = sem(raw_2*1000, nan_policy='omit')

    f, ax = plt.subplots()


    #Paula
    #plot bargraph
    ax.bar(1,mean_1,capsize=0.5, yerr=mean_sem_1, label='control', alpha=0.5, width=0.08, color='grey')
    ax.bar(1.1,mean_2, yerr=mean_sem_2, capsize=0.5, label='0,1µM Melatonin', alpha=0.7, width=0.08, color='orange')

    x_1 = [1] * raw_1.shape[0]
    x_2 = [1.1] * raw_2.shape[0]

    dots = ax.scatter(x_1, raw_1*1000, color = 'grey',  alpha = 0.4)
    jitter_dots(dots)
    dots = ax.scatter(x_2, raw_2*1000, color = 'orange', alpha = 0.4)
    jitter_dots(dots)

    ax.set_xlabel('Group')
    ax.set_ylabel('Mean time to first bout [ms]')
    ax.set_title('Time to first bout')
    #ax.set_ylim(0,1000)
    ax.legend()
    ax.grid(False)

    sns.despine(ax=ax, top=True, right=True)

    f.savefig(save_dir / f'fig_first_mean.pdf')

    plt.close(f)

    return 0

if __name__ == '__main__':

    experiment = 'd7_12_22_2021_0,1µM_Melatonin'
    #for experiment in ['d7_07_07_2021', 'd7_07_13_2021','d7_07_14_2021', 'd8_07_02_2021', 'd8_07_14_2021', 'd8_07_15_2021','d8_08_19_2021_recovery', 'd8_08_19_2021']:
    #main(experiment)
    plotHistogram(experiment)

    sys.exit(0)
