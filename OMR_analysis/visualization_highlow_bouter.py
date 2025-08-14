"""
run this code to plot the high & low bouters after having run all the other SD codes for all SD groups
"""
import numpy as np
import hdf5storage as hdf
import path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import sem

import pingouin as pg
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd
import scipy.stats as ss

def makeDf(data_1, data_2, stimuli, rename):
    '''
    if stimuli % 2 == 0:

        half = stimuli // 2

        print(data_1)
        data_1 = (data_1[:, :half] + data_1[:, half:]) / 2.0
        data_2 = (data_2[:, :half] + data_2[:, half:]) / 2.0

    else:

        half = (stimuli - 1) // 2

        data_1 = (data_1[:, :half] + data_1[:, half:-1]) / 2.0
        data_2 = (data_2[:, :half] + data_2[:, half:-1]) / 2.0
    '''
    stacked = np.concatenate((data_1, data_2), axis=0)
    col_names = [str(i) for i in range(stacked.shape[1])]

    df_temp = pd.DataFrame(stacked, columns=col_names)
    df_temp['Group'] = np.concatenate(([1] * data_1.shape[0], [2] * data_2.shape[0]))
    df_temp['ID'] = np.concatenate((np.arange(data_1.shape[0]), np.arange(data_2.shape[0])))

    df_data = df_temp.melt(id_vars=['Group', 'ID'], value_vars=col_names,
                           var_name='Stimulus', value_name=rename)

    # Paula: add new_ID column to identify fish based on stimulus value (indicating that the same fish has values for all stimuli)
    #fish_no =
    fish_list = range(0, 48)
    print(len(fish_list))
    print(len(df_data))
    n = 4
    repeated_list = []
    for _ in range(n):
        repeated_list.extend(fish_list)
    # print(repeated_list)
    df_data["New_ID"] = repeated_list
    # print(df_data)

    return df_data

def jitter_dots(dots, base_x = None, offsets = 1):
    offsets = dots.get_offsets()
    jittered_offsets = np.copy(offsets)
    if base_x is not None:
        jittered_offsets[:,0] = base_x
    # only jitter in the x-direction
    jittered_offsets[:, 0] += np.random.uniform(-0.05, 0.05, offsets.shape[0])
    dots.set_offsets(jittered_offsets)

def tfb(top, bottom, group, experiment):

    data_file = path.Path() / '..' / experiment
    tmp = hdf.loadmat(data_file/'data_all.mat')
    times = tmp[f'tfb_{group}']

    # Only 100% coherence for TFB

    stacked_top_times = np.stack((times[top][:,:,3], times[top][:,:,7]), axis=2)
    stacked_bottom_times = np.stack((times[bottom][:,:,3], times[bottom][:,:,7]), axis=2)

    top_times = np.nanmean(stacked_top_times, axis=(1,2)).reshape(-1)
    bottom_times = np.nanmean(stacked_bottom_times, axis=(1,2)).reshape(-1)

    sns.set_style('ticks')
    f, ax = plt.subplots()

    # Only 4 fish, so histogram is too sparse
    # Showing bar plot with the scattered individual values is better
    
    ax.set_ylabel('Reaction Time (s)')
    ax.grid(False)
    ax.set_ylim(0,max(bottom_times)+0.1)
    sns.despine(top=True, right=True)

    ax.bar([0.5, 1.5], [np.nanmean(top_times), np.nanmean(bottom_times)], \
           yerr=[np.nanstd(top_times), np.nanstd(bottom_times)], \
           tick_label=['Top', 'Bottom'], color=['royalblue', 'lightsteelblue'], alpha = 0.7)
    
    for i in range(top_times.shape[0]):
        #ax.scatter(0.5, top_times[i], c='black', alpha=0.5, s=10)
        #ax.scatter(1.5, bottom_times[i], c='black', alpha=0.5, s=10)

        dots = ax.scatter(0.5,top_times[i], color='royalblue', alpha=0.9)
        jitter_dots(dots)
        dots = ax.scatter(1.5,bottom_times[i], color='lightsteelblue', alpha=0.9)
        jitter_dots(dots)


    save_dir = path.Path() / '..' / experiment
    f.savefig(save_dir/f'tfb_bar_top_bottom.pdf')
    plt.close(f)

    # print means and SEM
    print('mean top', np.nanmean(top_times))
    print('STD top', np.nanstd(top_times))
    print('mean bottom', np.nanmean(bottom_times))
    print('STD bottom', np.nanstd(bottom_times))


    # Histogram with all individual values

    v_top, b_top = np.histogram(np.nanmean(stacked_top_times, axis=2), range=(0,2), bins=50)
    v_bottom, b_bottom = np.histogram(np.nanmean(stacked_bottom_times, axis=2), range=(0,2), bins=50)
    v_top = v_top / v_top.sum()
    v_bottom = v_bottom / v_bottom.sum()

    bins = (b_bottom[:-1] + b_bottom[1:]) / 2

    print(len(bins))

    f, ax = plt.subplots()

    ax.plot(bins, v_top, label='top', color='royalblue')
    ax.plot(bins, v_bottom, label='bottom', color='lightsteelblue')

    # from angles to calculate SEMs
    v_top_sem = sem(v_top, axis=0)
    v_bottom_sem = sem(v_bottom, axis=0)

    ax.fill_between(bins, v_top - v_top_sem, v_top + v_top_sem, color='royalblue', alpha=0.3)
    ax.fill_between(bins, v_bottom - v_bottom_sem, v_bottom+ v_bottom_sem,color='lightsteelblue', alpha=0.3)
    # end from angles

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Count')
    ax.set_xlim(0,2)
    #ax.set_ylim(0,0.11)
    ax.legend()
    ax.grid(False)
    sns.despine(top=True, right=True)

    save_dir = path.Path() / '..' / experiment
    f.savefig(save_dir/f'tfb_hist_top_bottom.pdf')
    plt.close(f)

    ## statistics
    # post hoc t-test for independent samples
    a = top_times
    b = bottom_times
    # for independent samples
    t_test = ss.ttest_ind(a, b)
    print(t_test)

    df_time = pd.DataFrame(data=t_test)

    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'time_stats.xlsx'
    df_time.to_excel(doc_name, index=False)

    return

def correct(top, bottom, group,experiment):

    data_file = path.Path() / '..' / experiment
    tmp = hdf.loadmat(data_file/'data_all.mat')
    performance = tmp[f'performance_{group}']

    top_perf_left = np.nanmean(performance[top][:,:,:4], axis=1)
    top_perf_right = np.nanmean(performance[top][:,:,4:], axis=1)
    bottom_perf_left = np.nanmean(performance[bottom][:,:,:4], axis=1)
    bottom_perf_right = np.nanmean(performance[bottom][:,:,4:], axis=1)

    top_perf = (top_perf_left + top_perf_right) / 2
    bottom_perf = (bottom_perf_left + bottom_perf_right) / 2

    print(top_perf)

    top_means = np.nanmean(top_perf, axis=0)
    bottom_means = np.nanmean(bottom_perf, axis=0)
    # top_std = np.nanstd(top_perf, axis=0)
    # bottom_std = np.nanstd(bottom_perf, axis=0)

    top_std = sem(top_perf, axis=0)
    bottom_std = sem(bottom_perf, axis=0)

    sns.set_style('ticks')
    f, ax = plt.subplots()

    offset = 0.03  # to distinguish groups more clearly

    ax.set_ylabel('Performance Score')
    ax.set_xlabel('Coherence')
    #ax.set_ylim(0.45,0.95)
    ax.grid(False)
    sns.despine(top=True, right=True)
    ax.set_xticks(np.arange(top_means.shape[0]), ['0', '25', '50', '100'])

    x_range = range(4)

    ax.errorbar(np.arange(top_means.shape[0]), top_means, yerr=top_std, color='royalblue', label='top')
    ax.errorbar(np.arange(top_means.shape[0]), bottom_means, yerr=bottom_std, color='lightsteelblue', label='bottom')

    for i in x_range:
        x_1 = np.full(top_perf.shape[0], i-offset)
        x_2 = np.full(bottom_perf.shape[0], i+ offset)

        dots = ax.scatter(x_1, top_perf[:, i], color='royalblue', alpha=0.3)
        jitter_dots(dots, base_x=x_1 - offset)
        dots = ax.scatter(x_2, bottom_perf[:, i], color='lightsteelblue', alpha=0.3)
        jitter_dots(dots, base_x=x_2 + offset)

    ax.legend()

    save_dir = path.Path() / '..' / experiment
    f.savefig(save_dir/f'correct_top_bottom.pdf')
    plt.close(f)

    ## statistics
    df_data = makeDf(top_perf,bottom_perf, stimuli, 'Correctness') #make data frame
    save_dir = path.Path() / '..' / experiment

    # mixed ANOVA (between and within for tests between different subjects who underwent several timepoints/conditions)
    results_mix = pg.mixed_anova(data=df_data, dv='Correctness', between='Group', within='Stimulus', subject='New_ID')
    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'correct_stats_mixed_anova.xlsx'
    results_mix.to_excel(doc_name, index=False)

    # Paula: posthoc t-test with scikit-posthoc
    # print("df_data",df_data.head())
    # for different stimuli
    df_stim_0 = df_data[df_data.Stimulus == "0"]  # replace number in "" for all stimuli (0-3)
    print("df_stim", df_stim_0)
    df_stim_1 = df_data[df_data.Stimulus == "1"]
    df_stim_2 = df_data[df_data.Stimulus == "2"]
    df_stim_3 = df_data[df_data.Stimulus == "3"]

    # posthoc tukey hsd
    posthoc_0 = pairwise_tukeyhsd(endog=df_stim_0['Correctness'], groups=df_stim_0['Group'], alpha=0.05)
    print("posthoc Correctness", posthoc_0)
    posthoc_1 = pairwise_tukeyhsd(endog=df_stim_1['Correctness'], groups=df_stim_1['Group'], alpha=0.05)
    posthoc_2 = pairwise_tukeyhsd(endog=df_stim_2['Correctness'], groups=df_stim_2['Group'], alpha=0.05)
    posthoc_3 = pairwise_tukeyhsd(endog=df_stim_3['Correctness'], groups=df_stim_3['Group'], alpha=0.05)

    # reuslts in dataframe
    df_correct_0 = pd.DataFrame(data=posthoc_0._results_table.data[1:], columns=posthoc_0._results_table.data[0])
    print(df_correct_0)
    df_correct_1 = pd.DataFrame(data=posthoc_1._results_table.data[1:], columns=posthoc_1._results_table.data[0])
    df_correct_2 = pd.DataFrame(data=posthoc_2._results_table.data[1:], columns=posthoc_2._results_table.data[0])
    df_correct_3 = pd.DataFrame(data=posthoc_3._results_table.data[1:], columns=posthoc_3._results_table.data[0])

    # Save in excel sheet
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(save_dir / 'correctness_posthoc_tukeyhsd.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df_correct_0.to_excel(writer, sheet_name='Stimuli 0%')
    df_correct_1.to_excel(writer, sheet_name='Stimuli 25%')
    df_correct_2.to_excel(writer, sheet_name='Stimuli 50%')
    df_correct_3.to_excel(writer, sheet_name='Stimuli 100%')
    # df3.to_excel(writer, sheet_name='Sheet3')

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    return


def angles(top, bottom, group, experiment):
    data_file = path.Path() / '..' / experiment
    tmp = hdf.loadmat(data_file/'data_all.mat')
    angles = tmp[f'angle_{group}']

    top_ang_left = np.nanmean(angles[top][:,:,:4,:], axis=1)
    top_ang_right = np.nanmean(angles[top][:,:,4:,:], axis=1)
    bottom_ang_left = np.nanmean(angles[bottom][:,:,:4,:], axis=1)
    bottom_ang_right = np.nanmean(angles[bottom][:,:,4:,:], axis=1)

    top_ang = (np.flip(top_ang_left, axis=2) + top_ang_right) / 2
    bottom_ang = (np.flip(bottom_ang_left, axis=2) + bottom_ang_right) / 2

    top_avg = np.nanmean(top_ang, axis=0)
    bottom_avg = np.nanmean(bottom_ang, axis=0)
    top_sem = sem(top_ang, axis=0)
    bottom_sem = sem(bottom_ang, axis=0)

    print(len(bottom_ang))

    angles = np.linspace(-180,180, angles.shape[3])
    sns.set_style('ticks')

    for stimulus in range(top_ang.shape[1]):

        f, ax = plt.subplots()

        ax.plot(angles, top_avg[stimulus], color='royalblue', label='top')
        ax.plot(angles, bottom_avg[stimulus], color='lightsteelblue', label='bottom')

        ax.fill_between(angles, top_avg[stimulus]-top_sem[stimulus],
                        top_avg[stimulus]+top_sem[stimulus],
                        color='royalblue', alpha=0.5)
        ax.fill_between(angles, bottom_avg[stimulus]-bottom_sem[stimulus],
                        bottom_avg[stimulus]+bottom_sem[stimulus],
                        color='lightsteelblue', alpha=0.5)

        ax.set_xlabel(f'$\\Delta$ Angle ($^\circ$)')
        ax.set_ylabel('Frequency (Hz)') # Verify and adjust data/units
        ax.set_title(f'Turn Distribution for {stimulus*25} %')
        ax.legend()
        ax.grid(False)
        ax.set_ylim(0,0.37)
        sns.despine(top=True, right=True)

        save_dir = path.Path() / '..' / experiment
        f.savefig(save_dir/f'angles_{stimulus}_top_bottom.pdf')
        plt.close(f)

    return


def getBouters(group, fraction, experiment):
    # Top and bottom bouters, each of size (1/fraction) * total fish
    data_file = path.Path() / '..'/experiment

    tmp = hdf.loadmat(data_file/'data_all.mat')
    rates = tmp[f'rates_{group}']

    rates_zero = np.stack((rates[:,:,0], rates[:,:,4]), axis=2).reshape(rates.shape[0], -1)
    rates_mean = np.nanmean(rates_zero, axis=1)
    sorted_rates = np.argsort(rates_mean)
    
    slice = int(sorted_rates.shape[0] // fraction)
    top = sorted_rates[-slice-1:]
    bottom = sorted_rates[:slice+1]

    return top, bottom

#test Paula
def boutFrequency(top, bottom, group,experiment, num_bins, stimuli):
    data_file = path.Path() / '..' / experiment
    tmp = hdf.loadmat(data_file / 'data_all.mat')
    rates = tmp[f'rates_{group}']

    save_dir = path.Path() / '..' / experiment

    top_rate_left = np.nanmean(rates[top][:, :, :4], axis=1)
    top_rate_right = np.nanmean(rates[top][:, :, 4:], axis=1)
    bottom_rate_left = np.nanmean(rates[bottom][:, :, :4], axis=1)
    bottom_rate_right = np.nanmean(rates[bottom][:, :, 4:], axis=1)

    top_rates = (top_rate_left + top_rate_right) / 2
    bottom_rates = (bottom_rate_left + bottom_rate_right) / 2
    print(top_rates)

    top_means = np.nanmean(top_rates, axis=0)
    bottom_means = np.nanmean(bottom_rates, axis=0)
    # top_std = np.nanstd(top_perf, axis=0)
    # bottom_std = np.nanstd(bottom_perf, axis=0)

    top_std = sem(top_rates, axis=0)
    bottom_std = sem(bottom_rates, axis=0)

    ## by Paula
    ## total response bout rate as a lineplot
    x_range = range(4)
    f, ax = plt.subplots()

    offset = 0.0

    ax.errorbar(np.arange(top_means.shape[0]), top_means, yerr=top_std, color='royalblue', label='top')
    ax.errorbar(np.arange(top_means.shape[0]), bottom_means, yerr=bottom_std, color='lightsteelblue', label='bottom')

    for i in x_range:
        x_1 = np.full(top_rates.shape[0], i - offset)
        x_2 = np.full(bottom_rates.shape[0], i + offset)

        dots = ax.scatter(x_1, top_rates[:, i], color='royalblue', alpha=0.3)
        jitter_dots(dots, base_x=x_1 - offset)
        dots = ax.scatter(x_2, bottom_rates[:, i], color='lightsteelblue', alpha=0.3)
        jitter_dots(dots, base_x=x_2 + offset)

    ax.legend()
    ax.set_xlabel(f'Coherence')
    ax.set_ylabel('Bouts/min')

    ax.set_title(f'Total response to stimulus')
    ax.set_xticks(x_range)
    text = [str(x) for x in [0, 25, 50, 100]]
    ax.set_xticklabels(text)
    ax.legend()
    #ax.set_ylim(0, 120)
    ax.grid(False)

    sns.despine(top=True, right=True)

    f.savefig(save_dir / f'fig_total_response_grey_doubled_lineplot_scatter.pdf')
    plt.close(f)

    # print mean, SEM, number of animals
    #print("Mean control:", freq_1[:4])
    #print("SEM control", sem_1)
    #print("Mean SD:", freq_2[:4])
    #print("SEM SD", sem_2)

    ## statistics
    df_data = makeDf(top_rates, bottom_rates, stimuli, 'Boutrates')  # make data frame
    save_dir = path.Path() / '..' / experiment

    # mixed ANOVA (between and within for tests between different subjects who underwent several timepoints/conditions)
    results_mix = pg.mixed_anova(data=df_data, dv='Boutrates', between='Group', within='Stimulus', subject='New_ID')
    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'boutrates_stats_mixed_anova.xlsx'
    results_mix.to_excel(doc_name, index=False)

    # Paula: posthoc t-test with scikit-posthoc
    # print("df_data",df_data.head())
    # for different stimuli
    df_stim_0 = df_data[df_data.Stimulus == "0"]  # replace number in "" for all stimuli (0-3)
    print("df_stim", df_stim_0)
    df_stim_1 = df_data[df_data.Stimulus == "1"]
    df_stim_2 = df_data[df_data.Stimulus == "2"]
    df_stim_3 = df_data[df_data.Stimulus == "3"]

    # posthoc tukey hsd
    posthoc_0 = pairwise_tukeyhsd(endog=df_stim_0['Boutrates'], groups=df_stim_0['Group'], alpha=0.05)
    print("posthoc Boutrates", posthoc_0)
    posthoc_1 = pairwise_tukeyhsd(endog=df_stim_1['Boutrates'], groups=df_stim_1['Group'], alpha=0.05)
    posthoc_2 = pairwise_tukeyhsd(endog=df_stim_2['Boutrates'], groups=df_stim_2['Group'], alpha=0.05)
    posthoc_3 = pairwise_tukeyhsd(endog=df_stim_3['Boutrates'], groups=df_stim_3['Group'], alpha=0.05)

    # reuslts in dataframe
    df_correct_0 = pd.DataFrame(data=posthoc_0._results_table.data[1:], columns=posthoc_0._results_table.data[0])
    print(df_correct_0)
    df_correct_1 = pd.DataFrame(data=posthoc_1._results_table.data[1:], columns=posthoc_1._results_table.data[0])
    df_correct_2 = pd.DataFrame(data=posthoc_2._results_table.data[1:], columns=posthoc_2._results_table.data[0])
    df_correct_3 = pd.DataFrame(data=posthoc_3._results_table.data[1:], columns=posthoc_3._results_table.data[0])

    # Save in excel sheet
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(save_dir / 'boutRates_posthoc_tukeyhsd.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df_correct_0.to_excel(writer, sheet_name='Stimuli 0%')
    df_correct_1.to_excel(writer, sheet_name='Stimuli 25%')
    df_correct_2.to_excel(writer, sheet_name='Stimuli 50%')
    df_correct_3.to_excel(writer, sheet_name='Stimuli 100%')
    # df3.to_excel(writer, sheet_name='Sheet3')

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    return

if __name__ == '__main__':

    experiment = 'highlow_bouter'
    group = 'control'
    fraction = 4 #  (1/fraction) * total fish. % fish you want to see for top/bottom 4 --> 25%
    info_path = path.Path() / '..' / experiment

    stimuli = 8
    num_bins = 72

    top, bottom = getBouters(group, fraction, experiment) # Find indices of top, bottom fish
    #tfb(top, bottom, group,experiment)
    #correct(top, bottom, group,experiment)
    #angles(top, bottom, group,experiment)
    boutFrequency(top, bottom, group, experiment, num_bins, stimuli)
    # correlation(group) # Not implemented ... yet