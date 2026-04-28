#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pyscript.py
#
#  Copyright 2020 Kumaresh <kumaresh_krishnan@g.harvard.edu>
#
#  version 2.0
#

import os, sys
import numpy as np

import hdf5storage as hdf
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mode

import path
import pickle

#Gaussian Fitting
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu, ttest_ind, sem
from itertools import combinations



def headingAngle(raw_data, stimulus, num_bins):
    start = 'bouts_start_stimulus_%03d' % (stimulus)
    end = 'bouts_end_stimulus_%03d' % (stimulus)

    # Compute differences (convention uses start-end)
    angles = raw_data[start]['fish_accumulated_orientation'] - \
             raw_data[end]['fish_accumulated_orientation']

    if angles.size == 0:
        return np.array([np.nan] * num_bins)

    # Find bout timestamps and fish location at start
    timestamps = raw_data[start]['timestamp']
    pos_x = raw_data[start]['fish_position_x']
    pos_y = raw_data[start]['fish_position_y']

    # Filter angles based on stimulus time and distance from edge of dish
    # Normalize by exposure to stimulus and suitable scaling for numerical value
    scale = 1000

    lim1, lim2 = 5.0, 15.00
    locs = np.where((timestamps > lim1) & (timestamps < lim2) & (pos_x ** 2 + pos_y ** 2 < 0.81))
    normalizer = (lim2 - lim1) / scale

    angles = angles[locs]

    # Restrict range to -180,180 and compute frequencies with specified num_bins
    angles[angles > 180] -= 360
    freq, _ = np.histogram(angles, bins=num_bins, range=(-180, 180))

    return freq / normalizer / 1000 * 60  # This is for bouts per minute


def extractAngles(experiment, root, num_bins):
    info_path = path.Path() / '..' / experiment

    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days = info['days']
    fish = info['fish']
    trials = info['trials']
    total_fish = np.sum(fish)

    fish_ctr = 0
    stimuli = 8

    data = np.full((total_fish, trials, stimuli, num_bins), np.nan)

    for day_idx, day in enumerate(days):

        for f in range(fish[day_idx]):

            for t in range(trials):

                folder = root / f'{day}_fish{f + 1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp = open(folder, 'rb')
                raw_data = pickle.load(tmp)

                for stimulus in range(stimuli):
                    angles = headingAngle(raw_data, stimulus, num_bins)
                    data[fish_ctr, t, stimulus] = angles

                tmp.close()

            fish_ctr += 1

        print(day, fish_ctr, 'fish done')

    return data


def processAngles(experiment, data, num_bins):
    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    group_1 = info['control']
    group_2 = info['sleep']

    avg_data_1 = np.nanmean(data[group_1], axis=(0, 1))
    avg_data_2 = np.nanmean(data[group_2], axis=(0, 1))

    sem_tmp_1 = np.nanmean(data[group_1], axis=1)
    sem_tmp_2 = np.nanmean(data[group_2], axis=1)

    sem_data_1 = sem(sem_tmp_1, axis=0, nan_policy='omit')
    sem_data_2 = sem(sem_tmp_2, axis=0, nan_policy='omit')

    norm_data_1 = avg_data_1.sum(axis=1).reshape(-1, 1)
    norm_data_2 = avg_data_2.sum(axis=1).reshape(-1, 1)

    prob_data_1 = avg_data_1 / norm_data_1
    prob_data_2 = avg_data_2 / norm_data_2

    prob_sem_1 = sem_data_1 / norm_data_1
    prob_sem_2 = sem_data_2 / norm_data_2

    tmp = np.nansum(data[group_1], axis=3)
    freq_tmp_1 = np.nanmean(tmp, axis=1)
    freq_data_1 = np.nanmean(freq_tmp_1, axis=0)
    freq_sem_1 = sem(freq_tmp_1, axis=0, nan_policy='omit')
    freq_std_1 = np.nanstd(freq_tmp_1, axis=0)

    tmp = np.nansum(data[group_2], axis=3)
    freq_tmp_2 = np.nanmean(tmp, axis=1)
    freq_data_2 = np.nanmean(freq_tmp_2, axis=0)
    freq_sem_2 = sem(freq_tmp_2, axis=0, nan_policy='omit')
    freq_std_2 = np.nanstd(freq_tmp_2, axis=0)

    to_save = {}
    to_save['mean_1'] = avg_data_1
    to_save['sem_1'] = sem_data_1.data if np.ma.isMaskedArray(sem_data_1) else sem_data_1
    to_save['prob_1'] = prob_data_1
    to_save['prob_sem_1'] = prob_sem_1.data if np.ma.isMaskedArray(prob_sem_1) else prob_sem_1
    to_save['mean_2'] = avg_data_2
    to_save['sem_2'] = sem_data_2.data if np.ma.isMaskedArray(sem_data_2) else sem_data_2
    to_save['prob_2'] = prob_data_2
    to_save['prob_sem_2'] = prob_sem_2.data if np.ma.isMaskedArray(prob_sem_2) else prob_sem_2
    # added for SEM
    to_save['freq_1_raw'] = freq_tmp_1
    to_save['freq_1'] = freq_data_1
    to_save['freq_sem_1'] = freq_sem_1.data if np.ma.isMaskedArray(freq_sem_1) else freq_sem_1
    to_save['freq_std_1'] = freq_std_1
    to_save['freq_2_raw'] = freq_tmp_2
    to_save['freq_2'] = freq_data_2
    to_save['freq_sem_2'] = freq_sem_2.data if np.ma.isMaskedArray(freq_sem_2) else freq_sem_2
    to_save['freq_std_2'] = freq_std_2

    return to_save


def main(experiment, num_bins):
    # root = path.Path() / '..' / '..' / '..' / 'data_hanna_test_06_16_2021' # directory for data
    root = path.Path() / '..' / experiment

    data = extractAngles(experiment, root, num_bins)

    to_save = processAngles(experiment, data, num_bins)

    save_dir = path.Path() / '..' / experiment / f'data_{num_bins}'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column', store_python_metadata=True)

    return 0


## add jitter: define jitter function
# https://stackoverflow.com/questions/60229586/how-to-make-jitterplot-on-matplolib-python
def jitter_dots(dots, base_x=None, offsets=1):
    offsets = dots.get_offsets()
    jittered_offsets = np.copy(offsets)
    if base_x is not None:
        jittered_offsets[:, 0] = base_x
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.05, 0.05, offsets.shape[0])
    dots.set_offsets(jittered_offsets)


def plotHistogram(experiment, num_bins, prob=False):
    data_path = path.Path() / '..' / experiment / f'data_{num_bins}'
    tmp = hdf.loadmat(data_path)

    if prob:
        data_1 = tmp['prob_1']
        data_2 = tmp['prob_2']
        sem_1 = tmp['prob_sem_1']
        sem_2 = tmp['prob_sem_2']
        save_dir = path.Path() / '..' / experiment / f'stimulus_histograms_{num_bins}_prob'
        save_dir_db = path.Path() / '..' / experiment / f'doubled_stimulus_histograms_{num_bins}_prob'

    else:
        data_1 = tmp['mean_1']
        data_2 = tmp['mean_2']
        sem_1 = tmp['sem_1']
        sem_2 = tmp['sem_2']
        save_dir = path.Path() / '..' / experiment / f'stimulus_histograms_{num_bins}'
        save_dir_db = path.Path() / '..' / experiment / f'doubled_stimulus_histograms_{num_bins}'

    stimuli = data_1.shape[0]

    angles = np.linspace(-180, 180, num_bins)
    id_map = hdf.loadmat(path.Path() / '..' / experiment / 'ID_map.mat')

    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(save_dir_db, exist_ok=True)
    sns.set()

    for stimulus in range(stimuli):

        f, ax = plt.subplots()

        ax.plot(angles, data_1[stimulus], label='control', color='black')
        ax.plot(angles, data_2[stimulus], label='light pulses morning', color='green')

        plt.fill_between(angles, \
                         data_1[stimulus] - sem_1[stimulus], \
                         data_1[stimulus] + sem_1[stimulus], \
                         color='grey', alpha=0.5)

        plt.fill_between(angles, \
                         data_2[stimulus] - sem_2[stimulus], \
                         data_2[stimulus] + sem_2[stimulus], \
                         color='green', alpha=0.5)

        ax.set_xlabel(f'$\\Delta$ Angle (°)')
        if prob == 1:
            ax.set_ylabel(f'Probability')
        else:
            ax.set_ylabel(f'Frequency (bpm)')  # Changed to bouts per minute
        ax.set_title(f'{id_map[str(stimulus)][0]} Stimulus {id_map[str(stimulus)][1]} %, (Bin: 5$^\circ$)')
        ax.legend()
        # sns.color_palette("PuRd", as_cmap=True)
        ax.grid(False)
        sns.set_style('white')
        sns.set_style('ticks')
        sns.despine(top=True, right=True)

        f.savefig(save_dir / f'fig_{stimulus}_{id_map[str(stimulus)][0]}_{id_map[str(stimulus)][1]}_grey.pdf')
        plt.close(f)

    if stimuli % 2 == 0:  # even numbers

        half = stimuli // 2

        data_1 = (np.fliplr(data_1[:half]) + data_1[half:]) / 2.0
        data_2 = (np.fliplr(data_2[:half]) + data_2[half:]) / 2.0

        sem_1 = np.sqrt((np.fliplr(sem_1[:half]) ** 2 + sem_1[half:] ** 2) / 2.0)
        sem_2 = np.sqrt((np.fliplr(sem_2[:half]) ** 2 + sem_2[half:] ** 2) / 2.0)

    else:  # odd numbers

        half = (stimuli - 1) // 2

        data_1 = (np.fliplr(data_1[:half]) + data_1[half:-1]) / 2.0
        data_2 = (np.fliplr(data_2[:half]) + data_2[half:-1]) / 2.0

        sem_1 = np.sqrt((np.fliplr(sem_1[:half]) ** 2 + sem_1[half:-1] ** 2) / 2.0)
        sem_2 = np.sqrt((np.fliplr(sem_2[:half]) ** 2 + sem_2[half:-1] ** 2) / 2.0)

    for stimulus in range(half):

        f, ax = plt.subplots()

        ax.plot(angles, data_1[stimulus], label='control', color='black')
        ax.plot(angles, data_2[stimulus], label='light pulses morning', color='green')

        plt.fill_between(angles, \
                         data_1[stimulus] - sem_1[stimulus], \
                         data_1[stimulus] + sem_1[stimulus], \
                         color='grey', alpha=0.5)

        plt.fill_between(angles, \
                         data_2[stimulus] - sem_2[stimulus], \
                         data_2[stimulus] + sem_2[stimulus], \
                         color='green', alpha=0.5)

        ax.set_xlabel(f'$\\Delta$ Angle (°)')
        if prob == 1:
            ax.set_ylabel(f'Probability')
        else:
            ax.set_ylabel(f'Frequency (bpm)')  # Changed to bouts per minute
        ax.set_title(
            f'{id_map[str(half + stimulus)][0]} Stimulus {id_map[str(half + stimulus)][1]} %, (Bin: 5$^\circ$)')
        ax.legend()
        ax.grid(False)
        sns.set_style('white')
        sns.set_style('ticks')
        ax.set_ylim(0, 0.30)
        sns.despine(top=True, right=True)

        # print(data_1)

        f.savefig(save_dir_db / f'fig_{id_map[str(half + stimulus)][0]}_{id_map[str(half + stimulus)][1]}_grey.pdf')
        plt.close(f)

        ## get angle distribution peaks

        # 1. Replace values below 10 with NaN
        data_1_filtered = np.where(data_1 >= 0.02, data_1, np.nan)
        print(data_1_filtered)

        # 2. Function to calculate the mode while ignoring NaNs
        def nan_mode(data):
            # Remove NaN values
            data = data[~np.isnan(data)]
            if len(data) == 0:
                return np.nan
            # Return the mode value
            mode_result = mode(data, nan_policy='omit')
            # if mode_result.count[0] > 0:
            # return mode_result.mode[0]
            # else:
            # return np.nan

        # 3. Calculate the median and mode for each stimulus
        median_data_1 = np.nanmedian(data_1_filtered, axis=1)
        mode_data_1 = np.array([nan_mode(row) for row in data_1_filtered])

        # Output the results
        print("Median values:", median_data_1)
        print("Mode values:", mode_data_1)
        '''
        # Filter the data, replacing values below 10 with NaN
        data_1_filtered = np.where(data_1 >= 10, data_1, np.nan)

        # Function to calculate mode while ignoring NaNs
        def nan_mode(data):
            # Remove NaN values
            data = data[~np.isnan(data)]
            if len(data) == 0:
                return np.nan
            # Return the mode value
            mode_result = mode(data, nan_policy='omit')
            if mode_result.count[0] > 0:
                return mode_result.mode[0]
            else:
                return np.nan

        # Calculate the median along axis 1 (for each stimulus)
        median_data_1 = np.array(
            [np.nanmedian(row) if not np.all(np.isnan(row)) else np.nan for row in data_1_filtered])

        # Calculate the mode along axis 1 (for each stimulus)
        mode_data_1 = np.array([nan_mode(row) for row in data_1_filtered])

        print("Median:", median_data_1)
        print("Mode:", mode_data_1)

        print("median control", median_data_1)
        print("median SD", median_data_2)
        print("mode control", mode_data_1)
        print("mode SD", mode_data_2)
        '''

        # data_1_filtered = np.where(data1_100 >= 10, data1_100, np.nan)
        # print(data_1_filtered)
        '''
        data1_stim0= data_1[stimulus==0]
        data1_stim25 = data_1[stimulus==1]
        data1_stim50 = data_1[stimulus==2]
        data1_stim100 = data_1[stimulus==3]

        data2_stim0 = data_2[stimulus==0]
        data2_stim25 = data_2[stimulus == 1]
        data2_stim50 = data_2[stimulus == 2]
        data2_stim100 = data_2[stimulus == 3]

        angle_dis_0 = data1_stim100[data1_stim100 > 3]
        print(angle_dis_0)
         '''
    return 0

"""
Gaussian peak fitting, extraction, and statistical comparison.
Append this section after the existing plotHistogram() function in pyscript.py.
Depends on: processAngles() output saved via hdf5storage, scipy, numpy, matplotlib, seaborn, pandas.
"""


# ---------------------------------------------------------------------------
# Gaussian model definitions
# c = centre, a = amplitude, w = width (sigma)
# ---------------------------------------------------------------------------

def gaussian_single(x, ctr, amp, wid):
    return amp * np.exp(-((x - ctr) / wid) ** 2)


def gaussian_2(x, c1, a1, w1, c2, a2, w2):
    return gaussian_single(x, c1, a1, w1) + gaussian_single(x, c2, a2, w2)


def gaussian_3(x, c1, a1, w1, c2, a2, w2, c3, a3, w3):
    return (gaussian_single(x, c1, a1, w1)
            + gaussian_single(x, c2, a2, w2)
            + gaussian_single(x, c3, a3, w3))


# ---------------------------------------------------------------------------
# fitGaussians
# ---------------------------------------------------------------------------

def fitGaussians(angles, y_data, n_gaussians=2, p0=None):
    """
    Fit n_gaussians Gaussians to a 1-D histogram.

    p0 format for 2 Gaussians: [c1, a1, w1, c2, a2, w2]
    p0 format for 3 Gaussians: [c1, a1, w1, c2, a2, w2, c3, a3, w3]

    If p0 is None a heuristic default is used — but you should always
    pass explicit p0 values from your __main__ block via p0_per_stimulus.
    """
    # Only use internal heuristic if caller did not supply p0
    if p0 is None:
        amp_est = float(np.nanmax(y_data))
        if n_gaussians == 2:
            # forward peak + right (correct) peak
            p0 = [0, amp_est * 0.8, 20,
                  45, amp_est * 0.5, 25]
        else:
            # left (indirect) + forward + right (correct)
            p0 = [-40, amp_est * 0.3, 20,
                  0, amp_est * 0.8, 20,
                  40, amp_est * 0.3, 20]

    func = gaussian_2 if n_gaussians == 2 else gaussian_3

    # Bounds: centres in [-180, 180], amplitudes >= 0, widths >= 1
    n_params = 6 if n_gaussians == 2 else 9
    lower = [-180, 0, 1] * n_gaussians
    upper = [180, np.inf, 180] * n_gaussians

    try:
        popt, pcov = curve_fit(
            func, angles, y_data,
            p0=p0, maxfev=10000,
            bounds=(lower, upper),
            method='trf'
        )
        fit = func(angles, *popt)
        return popt, pcov, fit, True
    except (RuntimeError, ValueError) as e:
        print(f"    curve_fit failed: {e}")
        nan_params = np.full(len(p0), np.nan)
        return nan_params, np.full((len(p0), len(p0)), np.nan), np.full_like(y_data, np.nan), False


# ---------------------------------------------------------------------------
# extractPeaks
# ---------------------------------------------------------------------------

def extractPeaks(popt, n_gaussians=2):
    """
    Returns list of dicts [{'centre', 'amplitude', 'width'}, ...] sorted
    left to right by centre.
    """
    peaks = []
    for i in range(n_gaussians):
        peaks.append({
            'centre': popt[i * 3],
            'amplitude': popt[i * 3 + 1],
            'width': popt[i * 3 + 2],
        })
    peaks.sort(key=lambda p: p['centre'])
    return peaks


# ---------------------------------------------------------------------------
# mirrorRaw  — apply once before any fitting
# ---------------------------------------------------------------------------

def mirrorRaw(raw):
    """
    Mirror and average stimulus pairs on the raw per-trial per-fish data.

    Input:  raw shape (n_fish, n_trials, n_stimuli, n_bins)  e.g. (N, T, 8, 72)
    Output: mirrored shape (n_fish, n_trials, n_stimuli//2, n_bins)  e.g. (N, T, 4, 72)

    Matches exactly what plotHistogram() does with np.fliplr on group means:
      mirrored[s] = (flipped(first_half[s]) + second_half[s]) / 2
    The flip is along the bin axis (axis=3), which reverses the angle direction.
    """
    n_stim = raw.shape[2]
    half = n_stim // 2 if n_stim % 2 == 0 else (n_stim - 1) // 2
    mirrored = (np.flip(raw[:, :, :half, :], axis=3)
                + raw[:, :, half:half + half, :]) / 2.0
    return mirrored  # (n_fish, n_trials, half, n_bins)


# ---------------------------------------------------------------------------
# fit_per_fish
# ---------------------------------------------------------------------------

def fit_per_fish(raw_mirrored, angles, n_gaussians_per_stimulus, p0_per_stimulus, peak_index):

    """
    For each fish and stimulus:
      1. Fit Gaussians to each individual trial
      2. Extract the requested peak centre from each trial
      3. Average those peak centres over trials → one value per fish per stimulus

    Parameters
    ----------
    raw_mirrored    : (n_fish, n_trials, n_stimuli, n_bins) — already mirrored
    angles          : (n_bins,) bin centres
    n_gaussians     : 2 or 3
    p0_per_stimulus : list of p0 guesses, one per stimulus (already group-specific)
    peak_index      : int (0=leftmost, 1=middle, 2=rightmost) or 'all' (mean of all)

    Returns
    -------
    peak_matrix : (n_fish, n_stimuli) — NaN where fitting failed for all trials
    """
    n_fish, n_trials, n_stim, _ = raw_mirrored.shape
    peak_matrix = np.full((n_fish, n_stim), np.nan)

    for f in range(n_fish):
        for s in range(n_stim):
            n_g = n_gaussians_per_stimulus[s]  # stimulus-specific
            p0 = (p0_per_stimulus[s]
                  if (p0_per_stimulus is not None and s < len(p0_per_stimulus))
                  else None)
            trial_peaks = []

            for t in range(n_trials):
                y = raw_mirrored[f, t, s, :]
                if np.all(np.isnan(y)):
                    continue
                popt, _, _, success = fitGaussians(angles, y, n_g, p0)
                if not success:
                    continue
                peaks = extractPeaks(popt, n_g)
                if peak_index == 'all':
                    trial_peaks.append(np.mean([p['centre'] for p in peaks]))
                elif isinstance(peak_index, int) and peak_index < len(peaks):
                    trial_peaks.append(peaks[peak_index]['centre'])

            if len(trial_peaks) > 0:
                peak_matrix[f, s] = np.nanmean(trial_peaks)

    return peak_matrix  # (n_fish, n_stimuli)


# ---------------------------------------------------------------------------
# runGaussianAnalysis  — top-level entry point
# ---------------------------------------------------------------------------

def runGaussianAnalysis(experiment, num_bins, n_gaussians_per_stimulus=None,
                        group1_label='control', group2_label='sleep',
                        peak_index=1,
                        p0_group1=None, p0_group2=None,
                        prob=True):
    """
    Full pipeline:
      mirror raw data → fit group means → fit per trial per fish →
      average peaks per fish → statistics → plots

    Parameters
    ----------
    experiment    : experiment folder name
    num_bins      : number of histogram bins (must match main())
    n_gaussians   : 2 or 3
    group1_label  : label for control group
    group2_label  : label for sleep/treatment group
    peak_index    : which Gaussian peak to compare:
                    0 = leftmost (indirect turns)
                    1 = middle (forward, if 3 Gaussians)
                    2 = rightmost (correct turns)
                    'all' = mean of all peak centres
    p0_group1     : list of p0 guesses per stimulus for group 1
                    e.g. [[c1,a1,w1, c2,a2,w2, c3,a3,w3], ...]  (one list per stimulus)
    p0_group2     : same format for group 2 (can differ from p0_group1)
    prob          : if True, work on probability data (values 0–1);
                    if False, work on raw frequencies (bpm)
    """

    # Default: 3 Gaussians for every stimulus if not specified
    if n_gaussians_per_stimulus is None:
        n_gaussians_per_stimulus = [3] * 4

    if peak_index is all or (not isinstance(peak_index, (int, str))):
        raise ValueError(
            "peak_index must be an integer or the string 'all'. "
            "Did you write peak_index=all instead of peak_index='all'?"
        )

    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    root = path.Path() / '..' / experiment
    angles = np.linspace(-180, 180, num_bins)

    group_1 = info['control']
    group_2 = info['sleep']

    # ------------------------------------------------------------------
    # 1. Load raw data and mirror ONCE
    # ------------------------------------------------------------------
    print("Loading raw data...")
    raw_data = extractAngles(experiment, root, num_bins)
    # raw_data shape: (n_fish_total, n_trials, 8_stimuli, n_bins)

    print("Mirroring stimulus pairs...")
    raw_g1 = mirrorRaw(raw_data[group_1])  # (n_fish1, n_trials, 4, n_bins)
    raw_g2 = mirrorRaw(raw_data[group_2])  # (n_fish2, n_trials, 4, n_bins)

    n_stimuli = raw_g1.shape[2]
    print(f"  Stimuli after mirroring: {n_stimuli}")

    # ------------------------------------------------------------------
    # 2. Compute group means + SEM for overlay plots
    #    average trials first, then fish
    # ------------------------------------------------------------------
    fish_mean_g1 = np.nanmean(raw_g1, axis=1)  # (n_fish1, 4, n_bins)
    fish_mean_g2 = np.nanmean(raw_g2, axis=1)  # (n_fish2, 4, n_bins)

    mean_g1 = np.nanmean(fish_mean_g1, axis=0)  # (4, n_bins)
    mean_g2 = np.nanmean(fish_mean_g2, axis=0)
    sem_g1 = sem(fish_mean_g1, axis=0, nan_policy='omit')
    sem_g2 = sem(fish_mean_g2, axis=0, nan_policy='omit')

    # If prob=True, normalise so each stimulus sums to 1
    if prob:
        norm1 = mean_g1.sum(axis=1, keepdims=True)
        norm2 = mean_g2.sum(axis=1, keepdims=True)
        # Avoid divide-by-zero
        norm1[norm1 == 0] = np.nan
        norm2[norm2 == 0] = np.nan
        mean_g1 = mean_g1 / norm1
        mean_g2 = mean_g2 / norm2
        sem_g1 = sem_g1 / norm1
        sem_g2 = sem_g2 / norm2
        ylabel = 'Probability'
        # Also normalise per-fish data for per-fish fitting
        norm_pf1 = fish_mean_g1.sum(axis=2, keepdims=True)
        norm_pf2 = fish_mean_g2.sum(axis=2, keepdims=True)
        norm_pf1[norm_pf1 == 0] = np.nan
        norm_pf2[norm_pf2 == 0] = np.nan
        # Normalise the raw trial data per trial
        norm_raw_g1 = raw_g1.sum(axis=3, keepdims=True)
        norm_raw_g2 = raw_g2.sum(axis=3, keepdims=True)
        norm_raw_g1[norm_raw_g1 == 0] = np.nan
        norm_raw_g2[norm_raw_g2 == 0] = np.nan
        raw_g1_fit = raw_g1 / norm_raw_g1  # each trial sums to 1
        raw_g2_fit = raw_g2 / norm_raw_g2
    else:
        ylabel = 'Frequency (bpm)'
        raw_g1_fit = raw_g1
        raw_g2_fit = raw_g2

    # ------------------------------------------------------------------
    # 3. Fit group means → overlay plots (4 stimuli only)
    # ------------------------------------------------------------------
    print(f"\n[1] Fitting group means ({n_gaussians_per_stimulus} Gaussians)...")
    id_map = hdf.loadmat(path.Path() / '..' / experiment / 'ID_map.mat')
    save_dir = path.Path() / '..' / experiment / f'gaussian_fits_{num_bins}'
    os.makedirs(save_dir, exist_ok=True)

    # fitAllStimuli uses whichever p0 you pass — group-specific
    fit_results_1 = fitAllStimuli(mean_g1, angles, n_gaussians_per_stimulus, p0_group1, group1_label)
    fit_results_2 = fitAllStimuli(mean_g2, angles, n_gaussians_per_stimulus, p0_group2, group2_label)

    all_peak_info = []

    for s in range(n_stimuli):
        # id_map keys 4–7 correspond to the mirrored (second-half) stimuli
        stim_key = str(n_stimuli + s)
        stim_name = f"{id_map[stim_key][0]} {id_map[stim_key][1]}%"
        r1 = fit_results_1[s]
        r2 = fit_results_2[s]

        fig, ax = plt.subplots(figsize=(7, 4))

        ax.plot(angles, mean_g1[s], color='black', label=group1_label, linewidth=1.2)
        ax.fill_between(angles,
                        mean_g1[s] - sem_g1[s],
                        mean_g1[s] + sem_g1[s],
                        color='grey', alpha=0.35)
        ax.plot(angles, mean_g2[s], color='red', label=group2_label, linewidth=1.2)
        ax.fill_between(angles,
                        mean_g2[s] - sem_g2[s],
                        mean_g2[s] + sem_g2[s],
                        color='red', alpha=0.25)

        for r, col in [(r1, 'black'), (r2, 'red')]:
            if r['success']:
                ax.plot(angles, r['fit'], color=col, linestyle='--',
                        linewidth=1.5, label=f'fit {col}')
                for peak in r['peaks']:
                    ax.axvline(peak['centre'], color=col, linestyle=':',
                               alpha=0.6, linewidth=1.0)
                    ax.text(peak['centre'] + 2,
                            ax.get_ylim()[1] * 0.92,
                            f"{peak['centre']:.1f}°",
                            color=col, fontsize=7)

        ax.set_xlabel('Δ Angle (°)')
        ax.set_ylabel(ylabel)
        ax.set_title(f'{stim_name} — {n_gaussians_per_stimulus}G fit, Bin: 5°')
        ax.legend(fontsize=8)
        sns.despine(top=True, right=True)
        fig.savefig(save_dir / f'fig_{s}_{stim_name}_gaussfit.pdf', bbox_inches='tight')
        plt.close(fig)

        for r, grp in [(r1, group1_label), (r2, group2_label)]:
            if r['success']:
                for pi, peak in enumerate(r['peaks']):
                    all_peak_info.append({
                        'stimulus': s, 'stim_label': stim_name, 'group': grp,
                        'peak_index': pi, 'centre_deg': peak['centre'],
                        'amplitude': peak['amplitude'], 'width': peak['width'],
                    })

    peak_df = pd.DataFrame(all_peak_info)
    peak_df.to_csv(save_dir / 'peak_summary_means.csv', index=False)
    print(f"  Saved peak_summary_means.csv")

    # ------------------------------------------------------------------
    # 4. Per-fish fitting:
    #    mirror done → fit each trial → average peaks over trials per fish
    #    Group-specific p0 passed explicitly
    # ------------------------------------------------------------------
    print(f"\n[2] Per-fish fitting ({group1_label})...")
    peaks_g1 = fit_per_fish(raw_g1_fit, angles, n_gaussians_per_stimulus, p0_group1, peak_index)

    print(f"\n[2] Per-fish fitting ({group2_label})...")
    peaks_g2 = fit_per_fish(raw_g2_fit, angles, n_gaussians_per_stimulus, p0_group2, peak_index)

    # ------------------------------------------------------------------
    # 5. Statistics: Mann-Whitney U + Welch's t-test per stimulus
    # ------------------------------------------------------------------
    print("\n[3] Statistics...")
    rows = []
    for s in range(n_stimuli):
        stim_key = str(n_stimuli + s)
        stim_name = f"{id_map[stim_key][0]} {id_map[stim_key][1]}%"

        v1 = peaks_g1[:, s][~np.isnan(peaks_g1[:, s])]
        v2 = peaks_g2[:, s][~np.isnan(peaks_g2[:, s])]

        row = {'stimulus': s, 'stim_label': stim_name,
               'n1': len(v1), 'n2': len(v2)}

        if len(v1) < 2 or len(v2) < 2:
            row.update({f'mean_{group1_label}': np.nan, f'sem_{group1_label}': np.nan,
                        f'mean_{group2_label}': np.nan, f'sem_{group2_label}': np.nan,
                        'mwu_stat': np.nan, 'mwu_p': np.nan,
                        'ttest_stat': np.nan, 'ttest_p': np.nan})
            rows.append(row)
            continue

        mwu_stat, mwu_p = mannwhitneyu(v1, v2, alternative='two-sided')
        tstat, tp = ttest_ind(v1, v2, equal_var=False)
        sig = ('***' if mwu_p < 0.001 else
               '**' if mwu_p < 0.01 else
               '*' if mwu_p < 0.05 else 'ns')

        print(f"  {stim_name}: {group1_label} μ={np.mean(v1):.1f}°  "
              f"{group2_label} μ={np.mean(v2):.1f}°  "
              f"MWU p={mwu_p:.4f} {sig}  Welch p={tp:.4f}")

        row.update({
            f'mean_{group1_label}': np.mean(v1),
            f'sem_{group1_label}': sem(v1),
            f'mean_{group2_label}': np.mean(v2),
            f'sem_{group2_label}': sem(v2),
            'mwu_stat': mwu_stat, 'mwu_p': mwu_p,
            'ttest_stat': tstat, 'ttest_p': tp,
        })
        rows.append(row)

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(save_dir / 'peak_stats.csv', index=False)
    print(f"  Saved peak_stats.csv")

    # ------------------------------------------------------------------
    # 6. Comparison plot
    # ------------------------------------------------------------------
    print("\n[4] Saving peak comparison plot...")
    plotPeakComparison(
        {'group1': peaks_g1, 'group2': peaks_g2},
        stats_df, experiment, num_bins,
        group1_label, group2_label
    )

    return peak_df, stats_df, {'group1': peaks_g1, 'group2': peaks_g2}


# ---------------------------------------------------------------------------
# fitAllStimuli — fits group-mean histograms, used by runGaussianAnalysis
# ---------------------------------------------------------------------------

def fitAllStimuli(data_matrix, angles, n_gaussians_per_stimulus, p0_per_stimulus, group_label):
    """
    data_matrix : (n_stimuli, n_bins) — group mean
    Returns list of result dicts, one per stimulus.
    """
    n_stimuli = data_matrix.shape[0]
    results = []
    for s in range(n_stimuli):
        n_g = n_gaussians_per_stimulus[s]  # stimulus-specific
        p0 = (p0_per_stimulus[s]
              if (p0_per_stimulus is not None and s < len(p0_per_stimulus))
              else None)
        popt, pcov, fit, success = fitGaussians(angles, data_matrix[s], n_g, p0)
        peaks = extractPeaks(popt, n_g) if success else []
        results.append({'stimulus': s, 'popt': popt, 'pcov': pcov,
                        'fit': fit, 'success': success, 'peaks': peaks})
        print(f"  [{group_label}] stimulus {s} ({n_g}G): {'OK' if success else 'FAILED'}")
    return results


# ---------------------------------------------------------------------------
# plotPeakComparison
# ---------------------------------------------------------------------------

def plotPeakComparison(fish_peaks, stats_df, experiment, num_bins,
                       group1_label='control', group2_label='sleep'):
    save_dir = path.Path() / '..' / experiment / f'gaussian_fits_{num_bins}'
    os.makedirs(save_dir, exist_ok=True)

    peaks1 = fish_peaks['group1']  # (n_fish1, n_stimuli)
    peaks2 = fish_peaks['group2']  # (n_fish2, n_stimuli)
    n_stimuli = peaks1.shape[1]

    fig, axes = plt.subplots(1, n_stimuli, figsize=(3 * n_stimuli, 4), sharey=True)
    if n_stimuli == 1:
        axes = [axes]

    def sig_label(p):
        if np.isnan(p): return ''
        if p < 0.001:   return '***'
        if p < 0.01:    return '**'
        if p < 0.05:    return '*'
        return 'ns'

    for s, ax in enumerate(axes):
        v1 = peaks1[:, s]
        v2 = peaks2[:, s]

        jitter1 = np.random.uniform(-0.08, 0.08, len(v1))
        jitter2 = np.random.uniform(-0.08, 0.08, len(v2))
        ax.scatter(np.zeros(len(v1)) + jitter1, v1, color='black', alpha=0.5, s=20, zorder=3)
        ax.scatter(np.ones(len(v2)) + jitter2, v2, color='red', alpha=0.5, s=20, zorder=3)

        for xi, (v, col) in enumerate([(v1, 'black'), (v2, 'red')]):
            valid = v[~np.isnan(v)]
            if len(valid) > 0:
                ax.errorbar(xi, np.mean(valid), yerr=sem(valid),
                            fmt='o', color=col, markersize=7,
                            capsize=4, linewidth=1.8, zorder=4)

        row = stats_df[stats_df['stimulus'] == s]
        if len(row):
            p = row['mwu_p'].values[0]
            label = sig_label(p)
            all_vals = np.concatenate([v1[~np.isnan(v1)], v2[~np.isnan(v2)]])
            if len(all_vals):
                ymax = np.max(all_vals) * 1.15
                ax.plot([0, 0, 1, 1],
                        [ymax * 0.92, ymax, ymax, ymax * 0.92],
                        lw=1, color='black')
                ax.text(0.5, ymax * 1.02, label,
                        ha='center', va='bottom', fontsize=10)

        stim_label = row['stim_label'].values[0] if len(row) else f'Stim {s}'
        ax.set_xticks([0, 1])
        ax.set_xticklabels([group1_label, group2_label], fontsize=8)
        ax.set_title(stim_label, fontsize=9)
        ax.axhline(0, color='grey', linewidth=0.6, linestyle='--', alpha=0.5)
        sns.despine(ax=ax, top=True, right=True)

    axes[0].set_ylabel('Peak angle (°)')
    fig.suptitle(f'Peak positions: {group1_label} vs {group2_label}', fontsize=11)
    fig.tight_layout()
    fig.savefig(save_dir / 'peak_comparison.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved peak_comparison.pdf")


# =============================================================================
# __main__
# =============================================================================

if __name__ == '__main__':
    experiment = 'd8_07_08_2021'
    num_bins = 72
    n_gaussians_per_stimulus = [3, 2, 2, 2]

    # -- Original pipeline (unchanged) ---------------------------------------
    main(experiment, num_bins)
    plotHistogram(experiment, num_bins)
    plotHistogram(experiment, num_bins, True)

    # -- Gaussian analysis ---------------------------------------------------
    # 9 values per stimulus: [c1, a1, w1,   c2, a2, w2,   c3, a3, w3]
    #                         left/indirect  forward       right/correct
    # Adjust centres and amplitudes to match what you see in your plots.

    p0_control = [
        [-30, 0.03, 10, 0, 0.25, 20, 30, 0.03, 10],  # stimulus 0 (lowest coherence)
        [ 0, 0.25, 20, 25, 0.05, 10],  # stimulus 1
        [ 0, 0.20, 12, 30, 0.05, 15],  # stimulus 2
        [ 0, 0.15, 10, 50, 0.10, 10],  # stimulus 3 (highest coherence)
    ]

    p0_sleep = [
        [-30, 0.03, 10, 0, 0.25, 20, 40, 0.03, 10],  # stimulus 0
        [ 0, 0.25, 20, 50, 0.05, 10],  # stimulus 1
        [ 0, 0.20, 10, 50, 0.05, 10],  # stimulus 2
        [ 0, 0.15, 10, 60, 0.10, 10],  # stimulus 3
    ]

    peak_df, stats_df, fish_peaks = runGaussianAnalysis(
        experiment, num_bins,
        n_gaussians_per_stimulus = n_gaussians_per_stimulus,  # one value per mirrored stimulus
        group1_label='control',
        group2_label='sleep',
        peak_index='all',  # 0=left/indirect  1=forward  2=right/correct
        p0_group1=p0_control,
        p0_group2=p0_sleep,
        prob=True,  # work on probability (0–1), not raw bpm
    )

    sys.exit(0)
