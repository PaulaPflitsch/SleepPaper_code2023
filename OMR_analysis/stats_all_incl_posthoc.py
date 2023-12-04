#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pyscript.py
#  
#  Copyright 2019 Kumaresh <kumaresh_krishnan@g.harvard.edu>
#  
#  version 1.0

import numpy as np
import os, sys
import csv

import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns

import hdf5storage as hdf
import path

import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from xlsxwriter import Workbook

from dfply import *


def makeDf(data_1, data_2, stimuli, rename):

    if stimuli % 2 == 0:

        half = stimuli // 2

        data_1 = (data_1[:,:half] + data_1[:,half:]) / 2.0
        data_2 = (data_2[:,:half] + data_2[:,half:]) / 2.0
        
    else:

        half = (stimuli - 1) // 2

        data_1 = (data_1[:,:half] + data_1[:,half:-1]) / 2.0
        data_2 = (data_2[:,:half] + data_2[:,half:-1]) / 2.0


    stacked = np.concatenate((data_1, data_2), axis=0)
    col_names = [str(i) for i in range(stacked.shape[1])]

    df_temp = pd.DataFrame(stacked, columns=col_names)
    df_temp['Group'] = np.concatenate(([1]*data_1.shape[0], [2]*data_2.shape[0]))
    df_temp['ID'] = np.concatenate((np.arange(data_1.shape[0]), np.arange(data_2.shape[0])))

    df_data = df_temp.melt(id_vars=['Group', 'ID'], value_vars=col_names,
        var_name='Stimulus', value_name=rename)

    #Paula: add new_ID column to identify fish based on stimulus value (indicating that the same fish has values for all stimuli)
    fish_list = range(0,64)
    print(len(fish_list))
    print(len(df_data))
    n = 4
    repeated_list = []
    for _ in range(n):
        repeated_list.extend(fish_list)
    #print(repeated_list)
    df_data["New_ID"]=repeated_list
    #print(df_data)

    return df_data

def totalBout(dpath, stimuli):

    tmp = hdf.loadmat(dpath / 'data_72.mat')
    data_1 = tmp['freq_1_raw']
    s=np.ndim(data_1)
    print(s)
    data_2 = tmp['freq_2_raw']

    df_data = makeDf(data_1, data_2, stimuli, 'Frequency')

    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'df_data.xlsx'
    df_data.to_excel(doc_name, index=False)

    bout_loc = df_data.loc[df_data['Stimulus'] =='0']
    sleep = bout_loc.loc[bout_loc['Group']== 2]
    control = bout_loc.loc[bout_loc['Group']==1]
    stats_bout = ss.ttest_ind(sleep.Frequency, control.Frequency, equal_var=False)

    sns.boxplot(x='Stimulus', y='Frequency', hue='Group', data=df_data, palette='Set3')

    # repeated measures ANOVA only for within fish controls, e.g. before and after treatment
    #results = pg.rm_anova(dv='Frequency', within=['Stimulus','Group'], subject='ID', data=df_data, detailed=True)
    #save_dir = path.Path() / '..' / experiment
    #doc_name = save_dir / 'bouts_stats.xlsx'
    #results.to_excel(doc_name, index=False)


    #mixed ANOVA (between and within for tests between different subjects who underwent several timepoints/condiitons)
    results_mix = pg.mixed_anova(data=df_data, dv='Frequency', between='Group', within='Stimulus', subject='New_ID')
    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'bouts_stats_mixed_anova.xlsx'
    results_mix.to_excel(doc_name, index=False)


    ## post hoc t-test for paired samples
    #print(df_data)
    #a = df_data.query('Group==1')['Frequency']
    #b = df_data.query('Group==2')['Frequency']
    #t_test = ss.ttest_ind(a, b)
    #print(t_test)

    # Paula: posthoc t-test with scikit-posthoc
    # print("df_data",df_data.head())
    # for different stimuli
    df_stim_0 = df_data[df_data.Stimulus == "0"]  # replace number in "" for all stimuli (0-3)
    print("df_stim", df_stim_0)
    df_stim_1 = df_data[df_data.Stimulus == "1"]
    df_stim_2 = df_data[df_data.Stimulus == "2"]
    df_stim_3 = df_data[df_data.Stimulus == "3"]

    # posthoc tukey hsd
    posthoc_0 = pairwise_tukeyhsd(endog=df_stim_0['Frequency'], groups=df_stim_0['Group'], alpha=0.05)
    print("posthoc Frequency", posthoc_0)
    posthoc_1 = pairwise_tukeyhsd(endog=df_stim_1['Frequency'], groups=df_stim_1['Group'], alpha=0.05)
    posthoc_2 = pairwise_tukeyhsd(endog=df_stim_2['Frequency'], groups=df_stim_2['Group'], alpha=0.05)
    posthoc_3 = pairwise_tukeyhsd(endog=df_stim_3['Frequency'], groups=df_stim_3['Group'], alpha=0.05)

    # reuslts in dataframe
    df_freq_0 = pd.DataFrame(data=posthoc_0._results_table.data[1:], columns=posthoc_0._results_table.data[0])
    df_freq_1 = pd.DataFrame(data=posthoc_1._results_table.data[1:], columns=posthoc_1._results_table.data[0])
    df_freq_2 = pd.DataFrame(data=posthoc_2._results_table.data[1:], columns=posthoc_2._results_table.data[0])
    df_freq_3 = pd.DataFrame(data=posthoc_3._results_table.data[1:], columns=posthoc_3._results_table.data[0])

    # Save in excel sheet
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(save_dir / 'Bout_freq_posthoc_tukeyhsd.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df_freq_0.to_excel(writer, sheet_name='Stimuli 0%')
    df_freq_1.to_excel(writer, sheet_name='Stimuli 25%')
    df_freq_2.to_excel(writer, sheet_name='Stimuli 50%')
    df_freq_3.to_excel(writer, sheet_name='Stimuli 100%')

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()


    return 0

def correctness(dpath, stimuli):

    tmp = hdf.loadmat(dpath / 'data_correctness_72.mat')
    data_1 = tmp['correct_1_raw']
    data_2 = tmp['correct_2_raw']

    df_data = makeDf(data_1, data_2, stimuli, 'Correctness')
    save_dir = path.Path() / '..' / experiment

    sns.boxplot(x='Stimulus', y='Correctness', hue='Group', data=df_data, palette='Set3')
    #plt.show()

    # anova between stimuli, between groups and between stimuli * groups for within subjects
    #results = pg.rm_anova(dv='Correctness', within=['Stimulus','Group'], subject='ID', data=df_data, detailed=True)
    #doc_name = save_dir / 'correct_stats.xlsx'
    #results.to_excel(doc_name, index=False)

    # mixed ANOVA (between and within for tests between different subjects who underwent several timepoints/condiitons)
    results_mix = pg.mixed_anova(data=df_data, dv='Correctness', between='Group', within='Stimulus', subject='New_ID')
    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'correct_stats_mixed_anova.xlsx'
    results_mix.to_excel(doc_name, index=False)


    # # ad hoc t-test for paired samples between groups
    a = df_data.query('Group==1')['Correctness']
    # print("a",a)
    b = df_data.query('Group==2')['Correctness']
    # print("b",b)
    # t_test = ss.ttest_ind(a,b)
    # print(t_test)


    # Paula: posthoc t-test with scikit-posthoc
    #print("df_data",df_data.head())
    # for different stimuli
    df_stim_0 = df_data[df_data.Stimulus=="0"] #replace number in "" for all stimuli (0-3)
    print("df_stim",df_stim_0)
    df_stim_1 = df_data[df_data.Stimulus == "1"]
    df_stim_2 = df_data[df_data.Stimulus == "2"]
    df_stim_3 = df_data[df_data.Stimulus == "3"]


    #posthoc tukey hsd
    posthoc_0 = pairwise_tukeyhsd(endog=df_stim_0['Correctness'],groups=df_stim_0['Group'],alpha=0.05)
    print("posthoc Correctness", posthoc_0)
    posthoc_1 = pairwise_tukeyhsd(endog=df_stim_1['Correctness'],groups=df_stim_1['Group'],alpha=0.05)
    posthoc_2 = pairwise_tukeyhsd(endog=df_stim_2['Correctness'], groups=df_stim_2['Group'], alpha=0.05)
    posthoc_3 = pairwise_tukeyhsd(endog=df_stim_3['Correctness'], groups=df_stim_3['Group'], alpha=0.05)

    # reuslts in dataframe
    df_correct_0 = pd.DataFrame(data=posthoc_0._results_table.data[1:], columns=posthoc_0._results_table.data[0])
    print(df_correct_0)
    df_correct_1 = pd.DataFrame(data=posthoc_1._results_table.data[1:], columns=posthoc_1._results_table.data[0])
    df_correct_2 = pd.DataFrame(data=posthoc_2._results_table.data[1:], columns=posthoc_2._results_table.data[0])
    df_correct_3 = pd.DataFrame(data=posthoc_3._results_table.data[1:], columns=posthoc_3._results_table.data[0])

    #Save in excel sheet
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(save_dir/'correctness_posthoc_tukeyhsd.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df_correct_0.to_excel(writer, sheet_name='Stimuli 0%')
    df_correct_1.to_excel(writer, sheet_name='Stimuli 25%')
    df_correct_2.to_excel(writer, sheet_name='Stimuli 50%')
    df_correct_3.to_excel(writer, sheet_name='Stimuli 100%')
    #df3.to_excel(writer, sheet_name='Sheet3')

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    return 0

def rate24(dpath):

    tmp = hdf.loadmat(dpath / 'data_rate.mat')
    data_1 = tmp['freq_1_raw']
    data_2 = tmp['freq_2_raw']

    g1_d1 = data_1[:,0:14*30].mean(axis=1) # test significance between 3 parts
    g1_night = data_1[:,14*30:34*30].mean(axis=1) # These have dimension fish x time series
    g1_d2 = data_1[:,35*30:48*30].mean(axis=1) # 2 dimension because trials are concatenated

    g2_d1 = data_2[:,0:14*30].mean(axis=1) # test significance between 3 parts
    g2_night = data_2[:,14*30:34*30].mean(axis=1) # These have dimension fish x time series
    g2_d2 = data_2[:,35*30:48*30].mean(axis=1) # 2 dimension because trials are concatenated

    stat_d1 = ss.ttest_ind(g1_d1, g2_d1, equal_var=False) #ad hoc independent t-test
    stat_night = ss.ttest_ind(g1_night, g2_night, equal_var=False)
    stat_d2 = ss.ttest_ind(g1_d2, g2_d2, equal_var=False)

    print(stat_d1); print(stat_night); print(stat_d2)
    # How do we want to save this? Excel? Will look into it and add
    
    return 0


def time(dpath):
    tmp = hdf.loadmat(dpath / 'data_time.mat')
    data_1 = tmp['raw_value_1']
    data_2 = tmp['raw_value_2']


    data_1_mean = np.nanmean(data_1, axis=(1, 2))
    data_2_mean = np.nanmean(data_2, axis=(1, 2))
    #data_1_mean = np.nanmean(data_1)
    #data_2_mean = np.nanmean(data_2)



    # post hoc t-test for independent samples or related samples
    #print(df_data)
    a = data_1_mean
    b = data_2_mean
    # for independent samples
    #t_test = ss.ttest_ind(a, b)
    # for related samples
    t_test = ss.ttest_rel(a, b)
    print(t_test)


    df_time = pd.DataFrame(data=t_test)


    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'time_stats.xlsx'
    df_time.to_excel(doc_name, index=False)
    return 0


def IBI(dpath, stimuli):
    tmp = hdf.loadmat(dpath / 'bout_data_50.mat')
    data_1 = tmp['raw_ib_1']
    data_2 = tmp['raw_ib_2']

    #Paula try:
    data_1 = np.nanmean(data_1, axis=2) # axis 2 is stimulus?!
    data_2 = np.nanmean(data_2, axis =2)

    df_data = makeDf(data_1, data_2, stimuli, 'IBI')

    bout_loc = df_data.loc[df_data['Stimulus'] == '0']
    sleep = bout_loc.loc[bout_loc['Group'] == 2]
    control = bout_loc.loc[bout_loc['Group'] == 1]
    #stats_bout = ss.ttest_ind(sleep.Frequency, control.Frequency, equal_var=False)

    sns.boxplot(x='Stimulus', y='IBI', hue='Group', data=df_data, palette='Set3')

    # for paired samples
    #results = pg.rm_anova(dv='IBI', within=['Stimulus', 'Group'], subject='ID', data=df_data, detailed=True)
    #save_dir = path.Path() / '..' / experiment
    #doc_name = save_dir / 'IBI_stats.xlsx'
    #results.to_excel(doc_name, index=False)

    # mixed ANOVA (between and within for tests between different subjects who underwent several timepoints/condiitons)
    results_mix = pg.mixed_anova(data=df_data, dv='IBI', between='Group', within='Stimulus', subject='New_ID')
    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'IBI_stats_mixed_anova.xlsx'
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
    posthoc_0 = pairwise_tukeyhsd(endog=df_stim_0['IBI'], groups=df_stim_0['Group'], alpha=0.05)
    print("posthoc Frequency", posthoc_0)
    posthoc_1 = pairwise_tukeyhsd(endog=df_stim_1['IBI'], groups=df_stim_1['Group'], alpha=0.05)
    posthoc_2 = pairwise_tukeyhsd(endog=df_stim_2['IBI'], groups=df_stim_2['Group'], alpha=0.05)
    posthoc_3 = pairwise_tukeyhsd(endog=df_stim_3['IBI'], groups=df_stim_3['Group'], alpha=0.05)

    # results in dataframe
    df_freq_0 = pd.DataFrame(data=posthoc_0._results_table.data[1:], columns=posthoc_0._results_table.data[0])
    df_freq_1 = pd.DataFrame(data=posthoc_1._results_table.data[1:], columns=posthoc_1._results_table.data[0])
    df_freq_2 = pd.DataFrame(data=posthoc_2._results_table.data[1:], columns=posthoc_2._results_table.data[0])
    df_freq_3 = pd.DataFrame(data=posthoc_3._results_table.data[1:], columns=posthoc_3._results_table.data[0])

    # Save in excel sheet
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(save_dir / 'IBI_posthoc_tukeyhsd.xlsx', engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    df_freq_0.to_excel(writer, sheet_name='Stimuli 0%')
    df_freq_1.to_excel(writer, sheet_name='Stimuli 25%')
    df_freq_2.to_excel(writer, sheet_name='Stimuli 50%')
    df_freq_3.to_excel(writer, sheet_name='Stimuli 100%')

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    return 0


def exploration(dpath):
    tmp = hdf.loadmat(dpath / 'data_extra_fishwise.mat')
    data_1 = tmp['raw_dist_1']
    print(data_1)
    data_2 = tmp['raw_dist_2']

    data_1_mean = np.nanmean(data_1, axis = 1)
    print("mean",data_1_mean)
    data_2_mean = np.nanmean(data_2, axis = 1)


    # post hoc t-test for independent samples or related samples
    # print(df_data)
    a = data_1_mean
    b = data_2_mean
    # for independent samples
    t_test = ss.ttest_ind(a,b)
    # for related samples
    #t_test = ss.ttest_rel(a, b)
    print(t_test)

    df_time = pd.DataFrame(data=t_test)

    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'exploration_stats.xlsx'
    df_time.to_excel(doc_name, index=False)

    return 0

def wrong_turns(dpath):
    tmp = hdf.loadmat(dpath / 'data_extra_fishwise.mat')
    data_1 = tmp['raw_wrong_1']
    print(data_1)
    data_2 = tmp['raw_wrong_2']

    data_1_mean = np.nanmean(data_1, axis = 1)
    print("mean",data_1_mean)
    data_2_mean = np.nanmean(data_2, axis = 1)


    # post hoc t-test for independent samples or related samples
    # print(df_data)
    a = data_1_mean
    b = data_2_mean
    # for independent samples
    t_test = ss.ttest_ind(a,b)
    # for related samples
    #t_test = ss.ttest_rel(a, b)
    print(t_test)

    df_time = pd.DataFrame(data=t_test)

    save_dir = path.Path() / '..' / experiment
    doc_name = save_dir / 'wrong_turns_stats.xlsx'
    df_time.to_excel(doc_name, index=False)

    return 0




if __name__ == '__main__':

    experiment = 'd7_07_07_2021'
    stimuli = 8

    dpath = path.Path() / '..' / experiment

    #totalBout(dpath, stimuli)
    #correctness(dpath, stimuli)
    #time(dpath)
    #IBI(dpath, stimuli)
    #exploration(dpath)
    wrong_turns(dpath)

    sys.exit()
