#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
import hdf5storage as hdf
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem, shapiro, mannwhitneyu
import path
import pickle


def headingAngle(raw_data, stimulus, num_bins):
    start = 'bouts_start_stimulus_%03d' % (stimulus)
    end = 'bouts_end_stimulus_%03d' % (stimulus)

    angles = raw_data[start]['fish_accumulated_orientation'] - \
             raw_data[end]['fish_accumulated_orientation']

    if angles.size == 0:
        return np.array([np.nan] * num_bins)

    timestamps = raw_data[start]['timestamp']
    pos_x = raw_data[start]['fish_position_x']
    pos_y = raw_data[start]['fish_position_y']

    scale = 1000
    lim1, lim2 = 5.0, 15.00
    locs = np.where((timestamps > lim1) & (timestamps < lim2) & (pos_x ** 2 + pos_y ** 2 < 0.81))
    normalizer = (lim2 - lim1) / scale

    angles = angles[locs]
    angles[angles > 180] -= 360

    freq, _ = np.histogram(angles, bins=num_bins, range=(-180, 180))
    return freq / normalizer / 1000 * 60


def extractRawAnglesForStimulus(raw_data, stimulus, mirror=False):
    """
    Extract filtered, wrapped raw delta angles and their timestamps for a
    single stimulus presentation within one trial.
    Timestamps are returned as time since stimulus onset (lim1 subtracted),
    so the axis runs 0–10s.

    mirror=True  : negate angles (leftward stimuli)
    mirror=False : keep angles as-is (rightward stimuli)
    """
    start = 'bouts_start_stimulus_%03d' % stimulus
    end   = 'bouts_end_stimulus_%03d' % stimulus

    angles = raw_data[start]['fish_accumulated_orientation'] - \
             raw_data[end]['fish_accumulated_orientation']

    if angles.size == 0:
        return np.array([]), np.array([])

    timestamps = raw_data[start]['timestamp']
    pos_x      = raw_data[start]['fish_position_x']
    pos_y      = raw_data[start]['fish_position_y']

    lim1, lim2 = 5.0, 15.00
    locs = np.where(
        (timestamps > lim1) & (timestamps < lim2) &
        (pos_x ** 2 + pos_y ** 2 < 0.81)
    )

    angles     = angles[locs]
    timestamps = timestamps[locs]
    angles[angles > 180] -= 360

    if mirror:
        angles = -angles

    sort_idx   = np.argsort(timestamps)
    # Subtract lim1 so timestamps represent time since stimulus onset (0–10s)
    return angles[sort_idx], timestamps[sort_idx] - lim1


def extractCumulativeByStimulus(experiment, root):
    """
    For each fish, each trial, and each of the 4 coherence levels:

    1. Compute the cumsum separately for each stimulus in the pair
       (leftward stimulus with sign flipped, rightward stimulus as-is).
    2. Store each as its own trace with real timestamps.

    The doubling (averaging the two paired traces) happens later in the
    plotting functions, after per-fish trial averaging.

    Data structure:
        doubled_fish[fish_idx][trial_idx][stim_pair_idx] = {
            'cumsum_first' : 1D cumsum for the leftward stimulus (sign flipped)
            'ts_first'     : timestamps for the leftward stimulus
            'cumsum_second': 1D cumsum for the rightward stimulus
            'ts_second'    : timestamps for the rightward stimulus
        }

    stim_pair_idx 0..3 corresponds to coherence 0%, 25%, 50%, 100%.
    The first half of the 8 stimuli (0..3) are leftward, the second half
    (4..7) are rightward — adjust mirror= if your convention differs.
    """
    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days        = info['days']
    fish_counts = info['fish']
    n_trials    = info['trials']
    group_1     = info['control']
    group_2     = info['sleep']

    stimuli = 8
    half    = stimuli // 2  # 4 coherence levels

    fish_ctr     = 0
    doubled_fish = []

    for day_idx, day in enumerate(days):
        for f in range(fish_counts[day_idx]):

            fish_trials = []

            for t in range(n_trials):
                folder   = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp      = open(folder, 'rb')
                raw_data = pickle.load(tmp)

                trial_pairs = []
                for s in range(half):
                    # Leftward stimulus: mirror=True so correct turns → positive
                    a_left,  ts_left  = extractRawAnglesForStimulus(
                        raw_data, s, mirror=True
                    )
                    # Rightward stimulus: mirror=False
                    a_right, ts_right = extractRawAnglesForStimulus(
                        raw_data, half + s, mirror=False
                    )

                    # Cumsum separately for each — resets to 0 for every
                    # stimulus presentation (they are sequential, independent)
                    cumsum_left  = np.cumsum(a_left)  if a_left.size  > 0 else np.array([np.nan])
                    cumsum_right = np.cumsum(a_right) if a_right.size > 0 else np.array([np.nan])

                    ts_left  = ts_left  if a_left.size  > 0 else np.array([np.nan])
                    ts_right = ts_right if a_right.size > 0 else np.array([np.nan])

                    trial_pairs.append({
                        'cumsum_left' : cumsum_left,
                        'ts_left'     : ts_left,
                        'cumsum_right': cumsum_right,
                        'ts_right'    : ts_right,
                    })

                tmp.close()
                fish_trials.append(trial_pairs)

            doubled_fish.append(fish_trials)
            fish_ctr += 1

        print(day, fish_ctr, 'fish done')

    return doubled_fish, group_1, group_2, half, n_trials


def _is_empty(arr):
    return arr.size == 0 or (arr.size == 1 and np.isnan(arr[0]))


def averageTraces(traces_list, common_axis, use_time=False):
    """
    Interpolate a list of (cumsum, x_axis) pairs onto common_axis and
    average, returning (mean, sem, valid_mask).

    traces_list : list of (cumsum_1d, x_axis_1d) tuples
    common_axis : 1D array — bout indices or time grid
    """
    matrix = np.full((len(traces_list), len(common_axis)), np.nan)

    for row, (c, x) in enumerate(traces_list):
        if _is_empty(c):
            continue
        in_range = (common_axis >= x[0]) & (common_axis <= x[-1])
        matrix[row, in_range] = np.interp(common_axis[in_range], x, c)

    mean  = np.nanmean(matrix, axis=0)
    s     = sem(matrix, axis=0, nan_policy='omit')
    valid = np.sum(~np.isnan(matrix), axis=0) >= 2

    return mean, s, valid, matrix


def buildPerFishTraces(doubled_fish, fish_indices, stim_idx, n_trials,
                       common_axis, use_time=False):
    """
    For each fish:
        1. For each trial, average the left and right cumulative traces
           onto common_axis (doubling step).
        2. Average those doubled traces across trials.

    Returns matrix of shape (n_fish, len(common_axis)).
    """
    matrix = np.full((len(fish_indices), len(common_axis)), np.nan)

    for row, fi in enumerate(fish_indices):
        trial_matrix = np.full((n_trials, len(common_axis)), np.nan)

        for t in range(n_trials):
            entry = doubled_fish[fi][t][stim_idx]

            c_left,  x_left  = entry['cumsum_left'],  entry['ts_left']
            c_right, x_right = entry['cumsum_right'], entry['ts_right']

            left_ok  = not _is_empty(c_left)
            right_ok = not _is_empty(c_right)

            if not left_ok and not right_ok:
                continue

            # Interpolate each available trace onto common_axis
            # then average them (doubling: two independent measurements
            # of the same coherence, now both positive = correct direction)
            interp_traces = []

            if left_ok:
                x = x_left if use_time else np.arange(len(c_left))
                in_range = (common_axis >= x[0]) & (common_axis <= x[-1])
                tmp = np.full(len(common_axis), np.nan)
                tmp[in_range] = np.interp(common_axis[in_range], x, c_left)
                interp_traces.append(tmp)

            if right_ok:
                x = x_right if use_time else np.arange(len(c_right))
                in_range = (common_axis >= x[0]) & (common_axis <= x[-1])
                tmp = np.full(len(common_axis), np.nan)
                tmp[in_range] = np.interp(common_axis[in_range], x, c_right)
                interp_traces.append(tmp)

            # Average left and right traces (doubling)
            trial_matrix[t] = np.nanmean(np.array(interp_traces), axis=0)

        # Average across trials for this fish
        matrix[row] = np.nanmean(trial_matrix, axis=0)

    return matrix


def plotCumulativeDeltaAngleByStimulus(experiment, root, stimulus_labels=None):
    """
    2x2 grid, one subplot per coherence level (0%, 25%, 50%, 100%).
    X-axis: bout index (starting from first bout in the stimulus window).
    For each fish: cumsum computed separately per stimulus direction per trial,
    the two directions averaged (doubling), then averaged across trials.
    Finally averaged across fish with SEM shading.
    """
    doubled_fish, group_1, group_2, half, n_trials = \
        extractCumulativeByStimulus(experiment, root)

    if stimulus_labels is None:
        stimulus_labels = [f'Stimulus {i}' for i in range(half)]

    # Common bout-index axis: find max bouts across all fish/trials/stimuli
    max_bouts = 1
    for fi in range(len(doubled_fish)):
        for t in range(n_trials):
            for s in range(half):
                entry = doubled_fish[fi][t][s]
                for key in ('cumsum_left', 'cumsum_right'):
                    c = entry[key]
                    if not _is_empty(c):
                        max_bouts = max(max_bouts, len(c))

    common_idx = np.arange(max_bouts)

    sns.set_style('white')
    sns.set_style('ticks')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
    axes = axes.flatten()

    for s in range(half):
        ax = axes[s]

        mat_1 = buildPerFishTraces(
            doubled_fish, group_1, s, n_trials, common_idx, use_time=False
        )
        mat_2 = buildPerFishTraces(
            doubled_fish, group_2, s, n_trials, common_idx, use_time=False
        )

        mean_1 = np.nanmean(mat_1, axis=0)
        mean_2 = np.nanmean(mat_2, axis=0)
        sem_1  = sem(mat_1, axis=0, nan_policy='omit')
        sem_2  = sem(mat_2, axis=0, nan_policy='omit')

        valid_1 = np.sum(~np.isnan(mat_1), axis=0) >= 2
        valid_2 = np.sum(~np.isnan(mat_2), axis=0) >= 2

        ax.plot(common_idx[valid_1], mean_1[valid_1],
                color='black', label='Control')
        ax.fill_between(common_idx[valid_1],
                        (mean_1 - sem_1)[valid_1],
                        (mean_1 + sem_1)[valid_1],
                        color='grey', alpha=0.4)

        ax.plot(common_idx[valid_2], mean_2[valid_2],
                color='red', label='Sleep Disruption')
        ax.fill_between(common_idx[valid_2],
                        (mean_2 - sem_2)[valid_2],
                        (mean_2 + sem_2)[valid_2],
                        color='red', alpha=0.3)

        ax.set_ylim(-160, 400)
        ax.axhline(90, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.text(max_bouts * 0.01, 92, '90°', fontsize=7, color='black', alpha=0.7)

        for mean, valid, color, name in [
            (mean_1, valid_1, 'black', 'Control'),
            (mean_2, valid_2, 'red',   'Sleep'),
        ]:
            valid_mean = mean[valid]
            valid_x    = common_idx[valid]
            crossings  = np.where(valid_mean >= 90)[0]
            if crossings.size > 0:
                cross_bout = valid_x[crossings[0]]
                ax.axvline(cross_bout, color=color,
                           linestyle='--', linewidth=0.8, alpha=0.7)
                ax.text(cross_bout + max_bouts * 0.01, 5,
                        f'b={cross_bout}', color=color, fontsize=7)
                print(f'Stimulus {s} | {name}: crosses 90° at bout {cross_bout}')

        ax.set_title(stimulus_labels[s])
        ax.set_xlabel('Bout Index')
        ax.set_ylabel('Cumulative Δ Angle (°)')
        ax.legend(fontsize=8)
        ax.grid(False)
        sns.despine(ax=ax, top=True, right=True)

    fig.suptitle('Cumulative Delta Angle by Stimulus (trial-averaged per fish)', y=1.01)
    fig.tight_layout()
    save_dir = path.Path() / '..' / experiment
    fig.savefig(save_dir / 'cumulative_delta_angle_by_stimulus_trial_avg.pdf',
                bbox_inches='tight')
    plt.close(fig)
    print('Saved cumulative_delta_angle_by_stimulus_trial_avg.pdf')


def plotCumulativeDeltaAngleOverTime(experiment, root, stimulus_labels=None):
    """
    2x2 grid, one subplot per coherence level.
    X-axis: real timestamps (ms since trial start, no subtraction).
    Same averaging hierarchy as the bout-index version, but interpolated
    onto a common time grid built from all real timestamps.
    """
    doubled_fish, group_1, group_2, half, n_trials = \
        extractCumulativeByStimulus(experiment, root)

    if stimulus_labels is None:
        stimulus_labels = [f'Stimulus {i}' for i in range(half)]

    # Build common time grid from union of all real timestamps
    all_ts = []
    for fi in range(len(doubled_fish)):
        for t in range(n_trials):
            for s in range(half):
                entry = doubled_fish[fi][t][s]
                for key in ('ts_left', 'ts_right'):
                    ts = entry[key]
                    if not _is_empty(ts):
                        all_ts.extend(ts.tolist())

    time_grid = np.array(sorted(set(all_ts)))

    # Time now runs 0-10s after stimulus onset sicne we earlier substracted 5 to account for the stimulus onset time
    lim1, lim2 = 0.0, 10.

    sns.set_style('white')
    sns.set_style('ticks')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=False, sharey=False)
    axes = axes.flatten()

    for s in range(half):
        ax = axes[s]

        mat_1 = buildPerFishTraces(
            doubled_fish, group_1, s, n_trials, time_grid, use_time=True
        )
        mat_2 = buildPerFishTraces(
            doubled_fish, group_2, s, n_trials, time_grid, use_time=True
        )

        mean_1 = np.nanmean(mat_1, axis=0)
        mean_2 = np.nanmean(mat_2, axis=0)
        sem_1  = sem(mat_1, axis=0, nan_policy='omit')
        sem_2  = sem(mat_2, axis=0, nan_policy='omit')

        valid_1 = np.sum(~np.isnan(mat_1), axis=0) >= 2
        valid_2 = np.sum(~np.isnan(mat_2), axis=0) >= 2

        ax.plot(time_grid[valid_1], mean_1[valid_1],
                color='black', label='Control')
        ax.fill_between(time_grid[valid_1],
                        (mean_1 - sem_1)[valid_1],
                        (mean_1 + sem_1)[valid_1],
                        color='grey', alpha=0.4)

        ax.plot(time_grid[valid_2], mean_2[valid_2],
                color='red', label='Sleep Disruption')
        ax.fill_between(time_grid[valid_2],
                        (mean_2 - sem_2)[valid_2],
                        (mean_2 + sem_2)[valid_2],
                        color='red', alpha=0.3)

        ax.set_ylim(-100, 400)
        ax.set_xlim(lim1, lim2)
        ax.axhline(90, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.text(lim1 + 0.05, 94, '90°', fontsize=7, color='black', alpha=0.7)

        for mean, valid, color, name in [
            (mean_1, valid_1, 'black', 'Control'),
            (mean_2, valid_2, 'red',   'Sleep'),
        ]:
            valid_mean = mean[valid]
            valid_time = time_grid[valid]
            crossings  = np.where(valid_mean >= 90)[0]
            if crossings.size > 0:
                cross_time = valid_time[crossings[0]]
                ax.axvline(cross_time, color=color,
                           linestyle='--', linewidth=0.8, alpha=0.7)
                ax.text(cross_time + (lim2 - lim1) * 0.01, 5,
                        f't={cross_time:.2f}s', color=color, fontsize=7)
                print(f'Stimulus {s} | {name}: crosses 90° at t={cross_time:.2f}s')

        ax.set_title(stimulus_labels[s])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Cumulative Δ Angle (°)')
        ax.legend(fontsize=8)
        ax.grid(False)
        sns.despine(ax=ax, top=True, right=True)

    fig.suptitle('Cumulative Delta Angle over Time (trial-averaged per fish)', y=1.01)
    fig.tight_layout()
    save_dir = path.Path() / '..' / experiment
    fig.savefig(save_dir / 'cumulative_delta_angle_over_time_avg.pdf',
                bbox_inches='tight')
    plt.close(fig)
    print('Saved cumulative_delta_angle_over_time_avg.pdf')


def get_crossing(cumsum, x_axis, threshold=90.0):
    """Return x value where cumsum first reaches threshold, else np.nan."""
    crossings = np.where(cumsum >= threshold)[0]
    return x_axis[crossings[0]] if crossings.size > 0 else np.nan


def statistics(experiment, root):
    """
    For each coherence level, test whether control and sleep groups differ
    in time (and bout index) to reach 90° cumulative delta angle.

    Per fish: crossing time/bout averaged across trials and across the two
    paired stimuli (left + right). Fish that never reach 90° are excluded.

    Normality tested with Shapiro-Wilk. Test selected automatically:
        - Both groups normal → independent t-test
        - Either group non-normal → Mann-Whitney U
    Bonferroni correction applied across the 4 coherence levels.
    """
    doubled_fish, group_1, group_2, half, n_trials = \
        extractCumulativeByStimulus(experiment, root)

    bonferroni_alpha = 0.05 / half

    print('\n=== Statistics: crossing 90° cumulative delta angle ===')
    print(f'Bonferroni-corrected alpha: {bonferroni_alpha:.4f}\n')

    for s in range(half):
        print(f'--- Coherence level {s} ---')

        for metric, use_time, unit in [
            ('Bout index', False, 'bouts'),
            ('Time',       True,  's'),
        ]:
            results = {}
            for fish_indices, name in [(group_1, 'Control'), (group_2, 'Sleep')]:
                per_fish = []
                for fi in fish_indices:
                    trial_crossings = []
                    for t in range(n_trials):
                        entry = doubled_fish[fi][t][s]
                        # Get crossing for left and right separately, average
                        pair_crossings = []
                        for c_key, ts_key in [
                            ('cumsum_left',  'ts_left'),
                            ('cumsum_right', 'ts_right'),
                        ]:
                            c  = entry[c_key]
                            ts = entry[ts_key]
                            if _is_empty(c):
                                continue
                            x_axis = ts if use_time else np.arange(len(c))
                            val = get_crossing(c, x_axis)
                            if not np.isnan(val):
                                pair_crossings.append(val)
                        if len(pair_crossings) > 0:
                            trial_crossings.append(np.mean(pair_crossings))
                    if len(trial_crossings) > 0:
                        per_fish.append(np.mean(trial_crossings))
                results[name] = np.array(per_fish)

            c1 = results['Control']
            c2 = results['Sleep']

            print(f'  [{metric}]')
            print(f'    Control : n={len(c1)}, mean={np.mean(c1):.3f} {unit}, '
                  f'sem={sem(c1):.3f}')
            print(f'    Sleep   : n={len(c2)}, mean={np.mean(c2):.3f} {unit}, '
                  f'sem={sem(c2):.3f}')

            sw_ok = True
            for arr, name in [(c1, 'Control'), (c2, 'Sleep')]:
                if len(arr) >= 3:
                    w, p_sw = shapiro(arr)
                    print(f'    Shapiro-Wilk {name}: W={w:.3f}, p={p_sw:.3f}')
                    if p_sw <= 0.05:
                        sw_ok = False
                else:
                    print(f'    Shapiro-Wilk {name}: n<3, skipped')
                    sw_ok = False

            if len(c1) >= 2 and len(c2) >= 2:
                if sw_ok:
                    from scipy.stats import ttest_ind
                    stat, p_val = ttest_ind(c1, c2)
                    test_name = 't-test'
                else:
                    stat, p_val = mannwhitneyu(c1, c2, alternative='two-sided')
                    test_name = 'Mann-Whitney U'
                sig = '***' if p_val < bonferroni_alpha else 'ns'
                print(f'    {test_name}: stat={stat:.3f}, p={p_val:.4f} {sig}')
            else:
                print(f'    Insufficient data for comparison')
        print()


def extractAngles(experiment, root, num_bins):
    info_path = path.Path() / '..' / experiment
    info = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
    days   = info['days']
    fish   = info['fish']
    trials = info['trials']
    total_fish = np.sum(fish)

    fish_ctr = 0
    stimuli  = 8
    data = np.full((total_fish, trials, stimuli, num_bins), np.nan)

    for day_idx, day in enumerate(days):
        for f in range(fish[day_idx]):
            for t in range(trials):
                folder   = root / f'{day}_fish{f+1:03d}' / 'raw_data' / f'trial{t:03d}.dat'
                tmp      = open(folder, 'rb')
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
    info    = np.load(info_path / 'expt_info.npy', allow_pickle=True).item()
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

    tmp        = np.nansum(data[group_1], axis=3)
    freq_tmp_1 = np.nanmean(tmp, axis=1)
    freq_data_1 = np.nanmean(freq_tmp_1, axis=0)
    freq_sem_1  = sem(freq_tmp_1, axis=0, nan_policy='omit')
    freq_std_1  = np.nanstd(freq_tmp_1, axis=0)

    tmp        = np.nansum(data[group_2], axis=3)
    freq_tmp_2 = np.nanmean(tmp, axis=1)
    freq_data_2 = np.nanmean(freq_tmp_2, axis=0)
    freq_sem_2  = sem(freq_tmp_2, axis=0, nan_policy='omit')
    freq_std_2  = np.nanstd(freq_tmp_2, axis=0)

    to_save = {
        'mean_1'    : avg_data_1,
        'sem_1'     : sem_data_1.data  if np.ma.isMaskedArray(sem_data_1)  else sem_data_1,
        'prob_1'    : prob_data_1,
        'prob_sem_1': prob_sem_1.data  if np.ma.isMaskedArray(prob_sem_1)  else prob_sem_1,
        'mean_2'    : avg_data_2,
        'sem_2'     : sem_data_2.data  if np.ma.isMaskedArray(sem_data_2)  else sem_data_2,
        'prob_2'    : prob_data_2,
        'prob_sem_2': prob_sem_2.data  if np.ma.isMaskedArray(prob_sem_2)  else prob_sem_2,
        'freq_1_raw': freq_tmp_1,
        'freq_1'    : freq_data_1,
        'freq_sem_1': freq_sem_1.data  if np.ma.isMaskedArray(freq_sem_1)  else freq_sem_1,
        'freq_std_1': freq_std_1,
        'freq_2_raw': freq_tmp_2,
        'freq_2'    : freq_data_2,
        'freq_sem_2': freq_sem_2.data  if np.ma.isMaskedArray(freq_sem_2)  else freq_sem_2,
        'freq_std_2': freq_std_2,
    }
    return to_save


def main(experiment, num_bins):
    root = path.Path() / '..' / experiment
    data = extractAngles(experiment, root, num_bins)
    to_save = processAngles(experiment, data, num_bins)
    save_dir = path.Path() / '..' / experiment / f'data_{num_bins}'
    hdf.savemat(save_dir, to_save, format='7.3', oned_as='column',
                store_python_metadata=True)
    return 0


if __name__ == '__main__':
    experiment = 'd7_07_07_2021'
    num_bins   = 72

    root = path.Path() / '..' / experiment

    main(experiment, num_bins)

    id_map = hdf.loadmat(path.Path() / '..' / experiment / 'ID_map.mat')
    labels = [f"{id_map[str(4+s)][0]} {id_map[str(4+s)][1]}%" for s in range(4)]

    plotCumulativeDeltaAngleByStimulus(experiment, root, stimulus_labels=labels)
    plotCumulativeDeltaAngleOverTime(experiment, root, stimulus_labels=labels)
    statistics(experiment, root)

    sys.exit(0)