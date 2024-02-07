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

# from matplotlib import cm
# from colorspacious import cspace_converter

def extraInfo(raw_data, stimulus):

    start = 'bouts_start_stimulus_%03d'%(stimulus)
    end = 'bouts_end_stimulus_%03d'%(stimulus)

    # Find bout timestamps and fish location at start
    timestamps = raw_data[start]['timestamp']
    pos_x = raw_data[end]['fish_position_x']
    pos_y = raw_data[end]['fish_position_y']
    angles = raw_data[start]['fish_accumulated_orientation'] - raw_data[end]['fish_accumulated_orientation']
    
    window = (timestamps > 5) & (timestamps < 15)
    if timestamps[window].size == 0:
        return np.nan, np.nan, np.nan, np.nan

    dist = np.sum(np.sqrt(np.diff(pos_x)**2 + np.diff(pos_y)**2)) # Difference from bout end

    react = timestamps[window][0] - 5

    if stimulus < 4:
        wrong = (angles[window] > 0).sum() + 0.5*(angles[window] == 0).sum()
        perf = angles[window][0] < 0 + 0.5*(angles[window][0] == 0)
    else:
        wrong = (angles[window] < 0).sum() + 0.5*(angles[window] == 0).sum()
        perf = angles[window][0] > 0 + 0.5*(angles[window][0] == 0)

    return react, wrong, dist, perf

def extractData(experiment, root):

    info_path = path.Path() / '..' / experiment

    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days = info['days']
    fish = info['fish']
    trials = info['trials']
    stimuli = 8 #info['stimuli'] #2
    print(stimuli)
    
    total_fish = np.sum(fish)
    
    fish_ctr = 0

    data_react = np.full((total_fish, trials, stimuli), np.nan)
    data_wrong = np.full((total_fish, trials, stimuli), np.nan)
    data_dist = np.full((total_fish, trials, stimuli), np.nan)
    data_perf = np.full((total_fish, trials, stimuli), np.nan)

    for day_idx, day in enumerate(days):

        for f in range(fish[day_idx]):

            for t in range(trials):

                for s_idx in range(stimuli): #all stimuli
                #for s_idx, s in enumerate([3,7]): # Only doing 100% coherence

                    folder = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                    tmp = open(folder, 'rb')
                    raw_data = pickle.load(tmp)
                        
                    react, wrong, dist, perf = extraInfo(raw_data, s_idx)
                    data_react[fish_ctr, t, s_idx] = react
                    data_wrong[fish_ctr, t, s_idx] = wrong
                    data_dist[fish_ctr, t, s_idx] = dist
                    data_perf[fish_ctr, t, s_idx] = perf

                tmp.close()

            fish_ctr += 1
                    
        print(day, fish_ctr, 'fish done')

    # Return a reshaped array with trials concatenated since it is 24 hr experiment
    
    return data_react, data_wrong, data_dist, data_perf

def processData(experiment, data_react, data_wrong, data_dist, data_perf):

    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    group_1 = info['control']
    group_2 = info['sleep']

    # sum across all stimuli
    #tot_dist_1 = np.nansum(data_dist[group_1], axis=2) # Sum across stimuli, average across trials
    #tot_dist_2 = np.nansum(data_dist[group_2], axis=2) # Sum across stimuli, average across trials

    # for stimulus 0% only

    tot_dist_1a= data_dist[group_1][:,:,0]
    tot_dist_1b = data_dist[group_1][:, :,4]

    tot_dist_1 = np.nansum(np.concatenate((tot_dist_1a[:,:,None], tot_dist_1b[:,:,None]), axis=2), axis=2)

    tot_dist_2a = data_dist[group_2][:, :, 0]
    tot_dist_2b = data_dist[group_2][:, :, 4]

    tot_dist_2 = np.nansum(np.concatenate((tot_dist_2a[:,:,None], tot_dist_2b[:,:,None]), axis=2), axis=2)

    #else:
        #tot_dist_1 = np.nansum(data_dist[group_1], axis=2)
        #tot_dist_2 = np.nansum(data_dist[group_2], axis=2)

    avg_dist_1 = np.nanmean(tot_dist_1, axis=1)
    avg_dist_2 = np.nanmean(tot_dist_2, axis=1)

    sem_dist_1 = sem(tot_dist_1, nan_policy='omit')
    sem_dist_2 = sem(tot_dist_2, nan_policy='omit')

    tot_wrong_1 = np.nansum(data_wrong[group_1], axis=2) # All wrong turns
    tot_wrong_2 = np.nansum(data_wrong[group_2], axis=2) # All wrong turns

    avg_wrong_1 = np.nanmean(tot_wrong_1, axis=1)
    avg_wrong_2 = np.nanmean(tot_wrong_2, axis=1)

    sem_wrong_1 = sem(tot_wrong_1, nan_policy='omit')
    sem_wrong_2 = sem(tot_wrong_2, nan_policy='omit')

    avg_react_1 = np.nanmean(data_react[group_1], axis=(1,2))
    avg_react_2 = np.nanmean(data_react[group_2], axis=(1,2))

    sem_react_1 = sem(data_react[group_1].reshape(len(group_1),-1), axis=1, nan_policy='omit')
    sem_react_2 = sem(data_react[group_2].reshape(len(group_2),-1), axis=1, nan_policy='omit')

    avg_perf_1 = np.nanmean(data_perf[group_1], axis=(1,2))
    avg_perf_2 = np.nanmean(data_perf[group_2], axis=(1,2))

    sem_perf_1 = sem(data_perf[group_1].reshape(len(group_1),-1), axis=1, nan_policy='omit')
    sem_perf_2 = sem(data_perf[group_2].reshape(len(group_2),-1), axis=1, nan_policy='omit')

    to_save = {}

    to_save['avg_dist_1'] = avg_dist_1
    to_save['sem_dist_1'] = sem_dist_1.data if np.ma.isMaskedArray(sem_dist_1) else sem_dist_1
    to_save['avg_wrong_1'] = avg_wrong_1
    to_save['sem_wrong_1'] = sem_wrong_1.data if np.ma.isMaskedArray(sem_wrong_1) else sem_wrong_1
    to_save['avg_react_1'] = avg_react_1
    to_save['sem_react_1'] = sem_react_1.data if np.ma.isMaskedArray(sem_react_1) else sem_react_1
    to_save['avg_perf_1'] = avg_perf_1
    to_save['sem_perf_1'] = sem_perf_1.data if np.ma.isMaskedArray(sem_perf_1) else sem_perf_1
    to_save['raw_dist_1'] = tot_dist_1
    to_save['raw_wrong_1'] = tot_wrong_1
    to_save['raw_react_1'] = data_react[group_1].reshape(len(group_1), -1)
    to_save['raw_perf_1'] = data_perf[group_1].reshape(len(group_1), -1)

    to_save['avg_dist_2'] = avg_dist_2
    to_save['sem_dist_2'] = sem_dist_2.data if np.ma.isMaskedArray(sem_dist_2) else sem_dist_2
    to_save['avg_wrong_2'] = avg_wrong_2
    to_save['sem_wrong_2'] = sem_wrong_2.data if np.ma.isMaskedArray(sem_wrong_2) else sem_wrong_2
    to_save['avg_react_2'] = avg_react_2
    to_save['sem_react_2'] = sem_react_2.data if np.ma.isMaskedArray(sem_react_2) else sem_react_2
    to_save['avg_perf_2'] = avg_perf_2
    to_save['sem_perf_2'] = sem_perf_2.data if np.ma.isMaskedArray(sem_perf_2) else sem_perf_2
    to_save['raw_dist_2'] = tot_dist_2
    to_save['raw_wrong_2'] = tot_wrong_2
    to_save['raw_react_2'] = data_react[group_2].reshape(len(group_2), -1)
    to_save['raw_perf_2'] = data_perf[group_2].reshape(len(group_2), -1)
    
    return to_save

def main(experiment):

    #root = path.Path() / '..' / '..' / '..' / 'data_hanna_test_06_16_2021' # directory for data
    root = path.Path() / '..' / experiment
    
    data_react, data_wrong, data_dist, data_perf = extractData(experiment, root,)

    to_save = processData(experiment, data_react, data_wrong, data_dist, data_perf)

    save_dir = path.Path() / '..' / experiment / f'data_extra_fishwise'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column', store_python_metadata=True)

    return 0

def jitter_dots(dots):
    offsets = dots.get_offsets()
    jittered_offsets = offsets
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.005, 0.005, offsets.shape[0])
    dots.set_offsets(jittered_offsets)

def plotWrong(experiment, prob=False):

    data_path = path.Path() / '..' / experiment / f'data_extra_fishwise'
    tmp = hdf.loadmat(data_path)

    data_1 = tmp['avg_wrong_1']
    data_2 = tmp['avg_wrong_2']
    sem_1 = tmp['sem_wrong_1']
    sem_2 = tmp['sem_wrong_2']

    save_dir = path.Path() / '..' / experiment

    sns.set_style('white')
    sns.set_style('ticks')    

    f, ax = plt.subplots()

    ax.bar(1, np.nanmean(data_1), \
        yerr=sem(data_1, nan_policy='omit'), width=0.08, \
        color='grey', capsize=0.5, alpha=0.5, label='control')

    ax.bar(1.1, np.nanmean(data_2), \
        yerr=sem(data_2, nan_policy='omit'), width=0.08, \
        color='orange', capsize=0.5, alpha=0.7, label='melatonin')

    x_1 = [1] * data_1.shape[0]
    x_2 = [1.1] * data_2.shape[0]

    dots = ax.scatter(x_1, data_1, color = 'grey',  alpha = 0.4)
    jitter_dots(dots)
    dots = ax.scatter(x_2, data_2, color = 'orange', alpha = 0.4)
    jitter_dots(dots)
    
    ax.set_xlabel(f'Group')
    ax.set_ylabel(f'Average Wrong turns')
    ax.set_title(f'#Wrong turns for each group')
    ax.legend()
    
    sns.color_palette("PuRd", as_cmap=True)
    ax.grid(False)

    sns.despine(top=True, right=True)

    f.savefig(save_dir / f'fig_wrong_turns_grey.pdf')
    plt.close(f)

    return 0

def plotDist(experiment, prob=False):

    data_path = path.Path() / '..' / experiment / f'data_extra_fishwise'
    tmp = hdf.loadmat(data_path)

    data_1 = tmp['avg_dist_1']
    data_2 = tmp['avg_dist_2']
    sem_1 = tmp['sem_dist_1']
    sem_2 = tmp['sem_dist_2']

    save_dir = path.Path() / '..' / experiment

    sns.set_style('white')
    sns.set_style('ticks')    

    f, ax = plt.subplots()

    ax.bar(1, np.nanmean(data_1), \
        yerr=sem(data_1, nan_policy='omit'), width=0.08, \
        color='grey', capsize=0.5, alpha=0.5, label='control')

    ax.bar(1.1, np.nanmean(data_2), \
        yerr=sem(data_2, nan_policy='omit'), width=0.08, \
        color='orange', capsize=0.5, alpha=0.7, label='melatonin')

    x_1 = [1] * data_1.shape[0]
    x_2 = [1.1] * data_2.shape[0]

    dots = ax.scatter(x_1, data_1, color = 'grey',  alpha = 0.4)
    jitter_dots(dots)
    dots = ax.scatter(x_2, data_2, color = 'orange', alpha = 0.4)
    jitter_dots(dots)
    
    ax.set_xlabel(f'Group')
    ax.set_ylabel(f'Average Distance')
    ax.set_title(f'Average distance covered by each group')
    ax.legend()
    
    sns.color_palette("PuRd", as_cmap=True)
    ax.grid(False)

    sns.despine(top=True, right=True)

    f.savefig(save_dir / f'fig_tot_dist_grey.pdf')
    plt.close(f)

    return 0

def plotReact(experiment, prob=False):

    data_path = path.Path() / '..' / experiment / f'data_extra_fishwise'
    tmp = hdf.loadmat(data_path)

    data_1 = tmp['avg_react_1']
    data_2 = tmp['avg_react_2']
    sem_1 = tmp['sem_react_1']
    sem_2 = tmp['sem_react_2']

    perf_1 = tmp['avg_perf_1']
    perf_2 = tmp['avg_perf_2']
    psem_1 = tmp['sem_perf_1']
    psem_2 = tmp['sem_perf_2']

    save_dir = path.Path() / '..' / experiment

    sns.set_style('white')
    sns.set_style('ticks')    

    f, ax = plt.subplots()

    ax.errorbar(data_1.ravel(), perf_1.ravel(), linestyle='', xerr=sem_1, yerr=psem_1, label='control', color='black', elinewidth = 0.2,marker='o', markersize = 3)
    ax.errorbar(data_2.ravel(), perf_2.ravel(),linestyle='', xerr=sem_2, yerr=psem_2, label='melatonin', color='orange', elinewidth = 0.2, marker='o' ,markersize =3 )
    
    ax.set_xlabel(f'Reaction time')
    ax.set_ylabel(f'Performance')
    ax.set_title(f'Reaction vs performance for first bout')
    ax.legend()
    #ax.set_xlim([0,3]) # Adjust
    ax.set_ylim([0,1.1])
    
    sns.color_palette("PuRd", as_cmap=True)
    ax.grid(False)

    sns.despine(top=True, right=True)

    f.savefig(save_dir / f'fig_reaction_grey.pdf')
    plt.close(f)

    return 0

if __name__ == '__main__':

    experiment = 'd7_12_22_2021_0,1µM_Melatonin'

    main(experiment)
    plotWrong(experiment)
    plotDist(experiment)
    plotReact(experiment)

    sys.exit(0)
