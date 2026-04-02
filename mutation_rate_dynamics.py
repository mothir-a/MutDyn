#--------------------------------------------------
# MRDAP: Mutation Rate Dynamics Analysis Program
# v4.0
# 20230120~
# Motohiro Akashi
#--------------------------------------------------

#--------------------------------------------------
# Import tools
#--------------------------------------------------

# set environment
import os
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import datetime
import gc

#--------------------------------------------------
# define initial dataframe
#--------------------------------------------------

# set initial df

def initial_population_df(initial_mut_rate_1, initial_mut_rate_2, initial_population_size):

    # set initial df
    df = pd.DataFrame(np.arange(initial_population_size * 11).reshape(initial_population_size, 11),
        columns=['gen', 'mut1', 'mut2', 'score_1', 'score_2', 'score_3', 'score_total', 'log10_mut_dif', 'mut_AVG', 'dead_alive', 'ancestor'])

    df['gen'] = 0
    df['mut1'] = float(initial_mut_rate_1)
    df['mut2'] = float(initial_mut_rate_2)
    df['score_1'] = 0
    df['score_2'] = 0
    df['score_3'] = 0
    df['score_total'] = df['score_1'] + df['score_2'] + df['score_3']
    df['log10_mut_dif'] = abs(np.log10(float(initial_mut_rate_1)) - np.log10(float(initial_mut_rate_2)))
    df['mut_AVG'] = np.log10(float((float(initial_mut_rate_1) + float(initial_mut_rate_2)) / 2 ))
    df['dead_alive'] = 1
    df['ancestor'] = df.index

    # set initial summary_df
    summary_df = df.describe().loc['mean':'mean']

    summary_df['cell_count'] = len(df)

    summary_df['survived_cell_count'] =  df['dead_alive'].sum()

    #return
    return df, summary_df

# END def initial_population_df-----------------------

#--------------------------------------------------
# define mutation rate dynamics function
#--------------------------------------------------

def mutation_rate_dynamics(df, summary_df, population_max, output_mid_data_interval, sharp_fall, sharp_fall_val, mut_dif_selection, mut_fix):
    #--------------------------------------------------
    # calculate by mut1
    #--------------------------------------------------

    # set df1
    df1 = df.copy()

    #(Note) If you write "df1=df", df will be changed when df1 is changed.

    df1['mut_rate'] = df1['mut1'].astype(float)

    #set ramd_choice column
    df1['rand_choice'] = ''

    #def random calculation
    def RandCal(row):
        row['rand_choice'] = np.random.choice(["a", "b", "c"], p=[row['mut_rate']/2, row['mut_rate']/2, 1 - row['mut_rate']])
        return row

    #calculate score_1
    df1 = df1.apply(RandCal, axis=1)
    df1.loc[df1['rand_choice'] == 'a','score_1'] -= 1
    df1.loc[df1['rand_choice'] == 'b','score_1'] += 1
    df1.loc[df1['rand_choice'] == 'c','score_1'] += 0

    #calculate score_2
    df1 = df1.apply(RandCal, axis=1)
    df1.loc[df1['rand_choice'] == 'a','score_2'] -= 1
    df1.loc[df1['rand_choice'] == 'b','score_2'] += 1
    df1.loc[df1['rand_choice'] == 'c','score_2'] += 0

    #calculate score_3
    df1 = df1.apply(RandCal, axis=1)
    df1.loc[df1['rand_choice'] == 'a','score_3'] -= 1
    df1.loc[df1['rand_choice'] == 'b','score_3'] += 1
    df1.loc[df1['rand_choice'] == 'c','score_3'] += 0

    #calculate mut2 (you must initially calculate mut2 mutation at mut1 stage)

    if mut_fix == 'on':
        df1['mut2'] = df1['mut2'].astype(float)
    elif mut_fix == 'off':
        df1 = df1.apply(RandCal, axis=1)
        df1.loc[df1['rand_choice'] == 'a','mut2'] *= 0.1
        df1.loc[df1['rand_choice'] == 'b','mut2'] /= 0.1
        df1.loc[df1['rand_choice'] == 'c','mut2'] *= 1
        df1['mut2'] = df1['mut2'].astype(float)

    #calculate mut1
    if mut_fix == 'on':
        df1['mut1'] = df1['mut1'].astype(float)
    elif mut_fix == 'off':
        df1 = df1.apply(RandCal, axis=1)
        df1.loc[df1['rand_choice'] == 'a','mut1'] *= 0.1
        df1.loc[df1['rand_choice'] == 'b','mut1'] /= 0.1
        df1.loc[df1['rand_choice'] == 'c','mut1'] *= 1
        df1['mut1'] = df1['mut1'].astype(float)


    #cleanup
    df1 = df1.drop(['rand_choice', 'mut_rate'], axis=1)

    #--------------------------------------------------
    # calculate by mut2
    #--------------------------------------------------

    # set df2
    df2 = df.copy()

    #(Note) If you write "df1=df", df will be changed when df1 is changed.
    df2['mut_rate'] = df2['mut2'].astype(float)

    #set ramd_choice column
    df2['rand_choice'] = ''

    #calculate score_1
    df2 = df2.apply(RandCal, axis=1)

    df2.loc[df2['rand_choice'] == 'a','score_1'] -= 1
    df2.loc[df2['rand_choice'] == 'b','score_1'] += 1
    df2.loc[df2['rand_choice'] == 'c','score_1'] += 0


    #calculate score_2
    df2 = df2.apply(RandCal, axis=1)

    df2.loc[df2['rand_choice'] == 'a','score_2'] -= 1
    df2.loc[df2['rand_choice'] == 'b','score_2'] += 1
    df2.loc[df2['rand_choice'] == 'c','score_2'] += 0

    #calculate score_3
    df2 = df2.apply(RandCal, axis=1)

    df2.loc[df2['rand_choice'] == 'a','score_3'] -= 1
    df2.loc[df2['rand_choice'] == 'b','score_3'] += 1
    df2.loc[df2['rand_choice'] == 'c','score_3'] += 0

    #calculate mut1 (you must initially calculate mut2 mutation at mut1 stage)

    if mut_fix == 'on':
        df2['mut1'] = df2['mut1'].astype(float)
    elif mut_fix == 'off':
        df2 = df2.apply(RandCal, axis=1)
        df2.loc[df2['rand_choice'] == 'a','mut1'] *= 0.1
        df2.loc[df2['rand_choice'] == 'b','mut1'] /= 0.1
        df2.loc[df2['rand_choice'] == 'c','mut1'] *= 1
        df2['mut1'] = df2['mut1'].astype(float)

    #calculate mut2
    if mut_fix == 'on':
        df2['mut2'] = df2['mut2'].astype(float)
    elif mut_fix == 'off':
        df2 = df2.apply(RandCal, axis=1)
        df2.loc[df2['rand_choice'] == 'a','mut2'] *= 0.1
        df2.loc[df2['rand_choice'] == 'b','mut2'] /= 0.1
        df2.loc[df2['rand_choice'] == 'c','mut2'] *= 1
        df2['mut2'] = df2['mut2'].astype(float)


    #cleanup
    df2 = df2.drop(['rand_choice', 'mut_rate'], axis=1)

    #--------------------------------------------------
    # concatenate df1 & df2
    #--------------------------------------------------

    # concatenate df1 & df2, output df_new
    df_new = pd.concat(
        [df1, df2],
        axis=0,
        ignore_index=True
    )

    #cleanup
    del df1, df2
    gc.collect()

    #--------------------------------------------------
    # write generation (+1)
    # Sharp Fall
    # sumup total score
    # calculate difference of 2 mutation rates
    # calculate average of 2 mutation rates (log10)
    # selection
    #--------------------------------------------------

    # write generation (+1)
    df_new['gen'] += 1

    # Sharp Fall
    if not sharp_fall == 'off' and not sharp_fall_val == 'off':
        if (df_new['gen'] == sharp_fall).any():
            df_new['score_1'] -= sharp_fall_val
            df_new['score_2'] -= sharp_fall_val
            df_new['score_3'] -= sharp_fall_val
        else:
            pass
    else:
        pass

    # sumup total score
    df_new['score_total'] = df_new['score_1'].astype(int) + df_new['score_2'].astype(int) + df_new['score_3'].astype(int)
    df_new['score_total'] = df_new['score_total'].astype(int)

    # calculate difference of 2 mutation rates
    df_new['log10_mut_dif'] = abs(np.log10(df_new['mut1'].astype(float)) - np.log10(df_new['mut2'].astype(float)))
    df_new['log10_mut_dif'] = df_new['log10_mut_dif'].astype(int)

    # calculate average of 2 mutation rates (log10)
    df_new['mut_AVG'] = np.log10((df_new['mut1'].astype(float) + df_new['mut2'].astype(float))/2)
    df_new['mut_AVG'] = df_new['mut_AVG'].astype(float)

    # selection
    # write down dead(0) or alive(2)
    df_new['selection_1'] = df_new['score_1']
    df_new['selection_2'] = df_new['score_2']
    df_new['selection_3'] = df_new['score_3']
    df_new['selection_4'] = abs(df_new['score_1'] - df_new['score_2'])
    df_new['selection_5'] = abs(df_new['score_2'] - df_new['score_3'])
    df_new['selection_6'] = df_new['mut1']
    df_new['selection_7'] = df_new['mut2']
    df_new['selection_8']  = df_new['log10_mut_dif']

    df_new['dead_alive'] = 1
    df_new.loc[df_new['selection_1'] < 0,'dead_alive'] = 0
    df_new.loc[df_new['selection_1'] > 100,'dead_alive'] = 0
    df_new.loc[df_new['selection_2'] < 0,'dead_alive'] = 0
    df_new.loc[df_new['selection_2'] > 100,'dead_alive'] = 0
    df_new.loc[df_new['selection_3'] < 0,'dead_alive'] = 0
    df_new.loc[df_new['selection_3'] > 100,'dead_alive'] = 0
    df_new.loc[df_new['selection_4'] > 2,'dead_alive'] = 0
    df_new.loc[df_new['selection_5'] > 2,'dead_alive'] = 0
    df_new.loc[df_new['selection_6'] > 1,'dead_alive'] = 0
    df_new.loc[df_new['selection_7'] > 1,'dead_alive'] = 0
    if not mut_dif_selection == 'off':
        df_new.loc[df_new['selection_8'] != int(mut_dif_selection),'dead_alive'] = 0


    #cleanup

    df_new = df_new.drop(['selection_1', 'selection_2', 'selection_3', 'selection_4', 'selection_5', 'selection_6', 'selection_7', 'selection_8'], axis=1)

    #--------------------------------------------------
    # sort & reindex df
    #--------------------------------------------------

    df_new = df_new.sort_values('score_total', ascending = False)

    df_new = df_new.reset_index(drop=True)

    #--------------------------------------------------
    # output df (before_selection)
    #--------------------------------------------------

    if df_new.describe().at['mean','gen'] == 1:
        #output csv
        gen_val = str(int(df_new.describe().at['mean','gen']))
        output_csv_adress = './mutation_rate_dynamics_gen' + gen_val + '.csv'
        df_new.to_csv(output_csv_adress, index=True, sep=',')


    if df_new.describe().at['mean','gen'] % output_mid_data_interval == 0:
        #output csv
        gen_val = str(int(df_new.describe().at['mean','gen']))
        output_csv_adress = './mutation_rate_dynamics_gen' + gen_val + '.csv'
        df_new.to_csv(output_csv_adress, index=True, sep=',')


    #--------------------------------------------------
    # sort & remove selected cells
    # cut off cells by population_max
    #--------------------------------------------------

    # sort & remove selected cells
    df_new_survive = df_new.sort_values('score_total', ascending = False)
    df_new_survive = df_new_survive[df_new_survive['dead_alive'] == 1]

    df_new_survive = df_new_survive.reset_index(drop=True)

    # cut off cells by population_max
    if len(df_new_survive) > population_max:
        df_new_survive = df_new_survive[0:population_max]

    #--------------------------------------------------
    # output df (survived)
    #--------------------------------------------------

    if df_new.describe().at['mean','gen'] == 1:
        #output csv
        gen_val = str(int(df_new.describe().at['mean','gen']))
        output_csv_adress = './mutation_rate_dynamics_gen' + gen_val + '_survived.csv'
        df_new_survive.to_csv(output_csv_adress, index=True, sep=',')


    if df_new.describe().at['mean','gen'] % output_mid_data_interval == 0:
        #output csv
        gen_val = str(int(df_new.describe().at['mean','gen']))
        output_csv_adress = './mutation_rate_dynamics_gen' + gen_val + '_survived.csv'
        df_new_survive.to_csv(output_csv_adress, index=True, sep=',')


    #--------------------------------------------------
    # summarizing
    #--------------------------------------------------

    # calculate mean values of each column.
    # Note: 'dead_alive' in summary_df indicates the survival rate.

    summary_df_prep = df_new_survive.describe().loc['mean':'mean']

    summary_df_prep['cell_count'] = len(df_new)

    summary_df_prep['survived_cell_count'] =  df_new_survive['dead_alive'].sum()


    summary_df = pd.concat(
        [summary_df, summary_df_prep],
        axis=0,
        ignore_index=True
    )

    #cleanup
    del summary_df_prep
    gc.collect()

    #--------------------------------------------------
    # overwrite df
    #--------------------------------------------------

    df = df_new_survive

    #cleanup
    del df_new, df_new_survive
    gc.collect()

    #--------------------------------------------------
    # return df
    #--------------------------------------------------

    return df, summary_df

# END def mutation_rate_dynamics-----------------------

#--------------------------------------------------
# define output summary graph
#--------------------------------------------------

def output_summary_fig(summary_df, fig_path, fig_dpi):

    # set dataset
    x = summary_df['gen'].astype(int)

    y1 = summary_df['score_total'].astype(float)
    y2 = summary_df['log10_mut_dif'].astype(float)
    y3 = summary_df['mut_AVG'].astype(float)

    # set plot area as a fig object
    fig = plt.figure(figsize = (6,10), facecolor='white')

    # set rcParams
    plt.rcParams.update({
        'font.family': 'Helvetica',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'legend.frameon': False,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'axes.linewidth': 0.5,
        'lines.linewidth': 1,
        'savefig.dpi': 450,
        'figure.figsize': (1.35, 3.35),  # 180 mm × 60 mm (Double column)
        'figure.dpi': 300
    })

    # set subplot area
    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    # set data to subplot areas
    ax1.plot(x, y1, color='#FEC10E', label='score_total')
    ax2.plot(x, y2, color='#FB5707', label='log10_mut_dif')
    ax3.plot(x, y3, color='#FD016E', label='mut_AVG')

    # add x label to subplot
    ax1.set_xlabel('Generation')
    ax2.set_xlabel('Generation')
    ax3.set_xlabel('Generation')

    # add y label to subplot
    ax1.set_ylabel('Total score')
    ax2.set_ylabel('$log_{10}$ average\nfidelity difference')
    ax3.set_ylabel('$log_{10}$ average\nmutation rate')

    # set grid
    ax1.grid(axis='both',linestyle='dotted', color='k')
    ax2.grid(axis='both',linestyle='dotted', color='k')
    ax3.grid(axis='both',linestyle='dotted', color='k')

    # spines invisible
    ax1.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # set legend
    #ax1.legend(loc = 'upper left')
    #ax2.legend(loc = 'upper left')
    #ax3.legend(loc = 'upper left')

    #save fig
    plt.savefig(fig_path, dpi=fig_dpi, bbox_inches='tight')

    #cleanup
    del x, y1, y2, y3, fig, ax1, ax2, ax3
    gc.collect()


# END def output_summary_fig------------------------------

#--------------------------------------------------
# define run_mutation_rate_dynamics
#--------------------------------------------------

def run_mutation_rate_dynamics(initial_mut_rate_1,
                               initial_mut_rate_2,
                               initial_population_size,
                               population_max,
                               output_mid_data_interval,
                               repetition_num,
                               sharp_fall,
                               sharp_fall_val,
                               mut_dif_selection,
                               mut_fix,
                               fig_dpi):

    #--------------------------------------------------
    # make & change directory
    #--------------------------------------------------

    # get datetime
    date_time = datetime.datetime.now()

    # datetime to string and set directory name
    date_time_str = date_time.strftime("%Y_%m%d_%H%M_%S")
    dir_name = 'test_' + date_time_str

    # make directory
    os.mkdir(dir_name)

    # change directory
    cd_path = './' + dir_name
    os.chdir(cd_path)

    #--------------------------------------------------
    # output initial setting
    #--------------------------------------------------

    # set initial_setting_df
    initial_setting_df = pd.DataFrame(columns=['dir_name',
                                               'cd_path',
                                               'initial_mut_rate_1',
                                               'initial_mut_rate_2',
                                               'initial_population_size',
                                               'population_max',
                                               'repetition_num',
                                               'sharp_fall',
                                               'sharp_fall_val',
                                               'mut_dif_selection',
                                               'mut_fix',
                                               'output_mid_data_interval',
                                               'summary_df_path',
                                               'summary_fig_path',
                                               'fig_dpi'],
                                      index = ['parameter'])

    # assign initial setting parameter
    initial_setting_df['dir_name'] = dir_name
    initial_setting_df['cd_path'] = cd_path
    initial_setting_df['initial_mut_rate_1'] = initial_mut_rate_1
    initial_setting_df['initial_mut_rate_2'] = initial_mut_rate_2
    initial_setting_df['initial_population_size'] = initial_population_size
    initial_setting_df['population_max'] = population_max
    initial_setting_df['repetition_num'] = repetition_num
    initial_setting_df['sharp_fall'] = sharp_fall
    initial_setting_df['sharp_fall_val'] = sharp_fall_val
    initial_setting_df['mut_dif_selection'] = mut_dif_selection
    initial_setting_df['mut_fix'] = mut_fix
    initial_setting_df['output_mid_data_interval'] = output_mid_data_interval
    initial_setting_df['summary_df_path'] = './mutation_rate_dynamics_summary_' + date_time_str + '.csv'
    initial_setting_df['summary_fig_path'] = './summary_df_fig_' + date_time_str + '.png'
    initial_setting_df['fig_dpi'] = fig_dpi

    # transpose initial_setting_df
    initial_setting_df = initial_setting_df.T

    # output initial_setting_df
    initial_setting_df_path = './initial_setting_' + date_time_str + '.csv'
    initial_setting_df.to_csv(initial_setting_df_path, index=True, sep=',')

    #--------------------------------------------------
    # set initial dataframe
    #--------------------------------------------------

    #initial setting for 'initial_population_df'
    initial_mut_rate_1 = initial_mut_rate_1
    initial_mut_rate_2 = initial_mut_rate_2
    initial_population_size = initial_population_size

    df, summary_df = initial_population_df(initial_mut_rate_1, initial_mut_rate_2, initial_population_size)

    #--------------------------------------------------
    # run mutation_rate_dynamics
    #--------------------------------------------------

    #initial setting for 'mutation_rate_dynamics' & number of calculation
    population_max = population_max
    output_mid_data_interval = output_mid_data_interval
    repetition_num = repetition_num
    sharp_fall = sharp_fall
    mut_dif_selection = mut_dif_selection
    mut_fix = mut_fix

    for value in range(repetition_num):
        df, summary_df = mutation_rate_dynamics(df, summary_df, population_max, output_mid_data_interval, sharp_fall, sharp_fall_val, mut_dif_selection, mut_fix)
        if df.empty == True:
            break

    #--------------------------------------------------
    # output summary_csv
    #--------------------------------------------------

    summary_df_path = './mutation_rate_dynamics_summary_' + date_time_str + '.csv'
    summary_df.to_csv(summary_df_path, index=True, sep=',')

    #--------------------------------------------------
    # output　summary graph
    #--------------------------------------------------

    # set fig path and resolution
    summary_fig_path = './summary_df_fig_' + date_time_str + '.png'
    fig_dpi = 200

    # output summary fig
    summary_df = summary_df
    output_summary_fig(summary_df, summary_fig_path, fig_dpi)

    #--------------------------------------------------
    # cleanup
    #--------------------------------------------------
    #cleanup
    del date_time, dir_name, cd_path, date_time_str
    del initial_setting_df
    del initial_mut_rate_1, initial_mut_rate_2, initial_population_size
    del population_max, output_mid_data_interval, repetition_num, sharp_fall, sharp_fall_val, mut_dif_selection, mut_fix
    del summary_df_path
    del summary_fig_path, fig_dpi
    del df, summary_df
    gc.collect()

    #--------------------------------------------------
    # change directory
    #--------------------------------------------------

    os.chdir('../')

# END def run_mutation_rate_dynamics----------------------


#--------------------------------------------------
#--------------------------------------------------
# test
#--------------------------------------------------
#--------------------------------------------------

# Input Initial Parameters

initial_mut_rate_1 = 0.1              # mut rate 1
initial_mut_rate_2 = 0.1              # mut rate 2
initial_population_size = 10        # initial population size
population_max = 50                 # under 5000 is recommended
output_mid_data_interval = 100      # output mid df per generation
repetition_num = 1000               # generation
sharp_fall = 'off'                    # generation in which sharp fall occur, 'off' or int
sharp_fall_val = 'off'                # sharp fall score, 'off' or int, lower than 99 is recommended
mut_dif_selection = 'off'             # selection by log10 difference of mut rates, 'off' or int
mut_fix = 'off'                       # fix mutation rate, 'on' or 'off'
fig_dpi = 200                         # summary fig dpi setting

# Execute run_mutation_rate_dynamics()
if __name__ == "__main__":

    run_mutation_rate_dynamics(initial_mut_rate_1,
                               initial_mut_rate_2,
                               initial_population_size,
                               population_max,
                               output_mid_data_interval,
                               repetition_num,
                               sharp_fall,
                               sharp_fall_val,
                               mut_dif_selection,
                               mut_fix,
                               fig_dpi)

#--------------------------------------------------
#--------------------------------------------------
# test END
#--------------------------------------------------
#--------------------------------------------------

###################################################
