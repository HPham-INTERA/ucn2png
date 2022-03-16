from funcs_ECF21_0006 import coc_prop
from funcs import read_ref, arr2shp, read_mas_files, generate_map1, func_compare_2ucn_files
import datetime as dt
import geopandas as gpd
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import pandas as pd
# !/opt/anaconda3.7/bin/python
# from mpl_toolkits.basemap import Basemap
# from pyproj import Proj, transform


# Run in cluster
# cd /workspace2/hpham/032_ECF_Hanford_21_0006/scripts
# ssh hpham@172.16.100.240
# Check mass error for CIE NFA (ECF-HANFORD-21-0006)
# 5/18/2021: Add processing *.mas
# 10/25/2021: This CIE MA script was modified in this version for CA Rev.1
# Last update: 12/09/2021


def var_unit(var):
    dic_units = {'cyan': '$\mu$g', 'cn': '$\mu$g',
                 'crvi': '$\mu$g', 'cr': '$\mu$g',
                 'i129': 'pCi', 'i-129': 'pCi',
                 'no3_': '$\mu$g', 'no3': '$\mu$g',
                 'sr90': 'pCi', 'sr-90': 'pCi',
                 'tc99': 'pCi', 'tc-99': 'pCi',
                 'tce': '$\mu$g',
                 'trit': 'pCi', 'h-3': 'pCi',
                 'utot': '$\mu$g', 'u-tot': '$\mu$g',
                 'u': '$\mu$g', 'ccl4': '$\mu$g'}
    # unit for var
    var_unit = dic_units[var]
    if var_unit == '$\mu$g':
        unit_conv_coef = 1e9  # from ug to kg
    elif var_unit == 'pCi':
        unit_conv_coef = 1e12
    return dic_units[var], unit_conv_coef


def create_new_dir(directory):
    # directory = os.path.dirname(file_path)
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
        print(f'Created a new directory {directory}\n')


def read_hss_files(path_to_files, var):
    '''
    hss files units: Time in days, mass_rate in
        - pCi /d if it is a Rad
        - ug if it is a chem
    '''
    os.chdir(path_to_files)
    # print(f'path_to_files@line53={path_to_files}\n')
    list_ifiles = glob("i*.dat")
    # print(f' Line 55, list_ifiles={list_ifiles}\n ')
    # unit, unit_conv_coef = var_unit(var) # CIE
    unit, unit_conv_coef = 'pCi', 1e12  # CA, radionuclides

    # get ijk for hss files
    dfMetaHss = pd.DataFrame(columns=['Name', 'I', 'J', 'K'])
    dfMetaHss['Name'] = list_ifiles
    for i, filename in enumerate(list_ifiles):
        row = int(filename.split('i')[1].split('j')[0])
        col = int(filename.split('i')[1].split('j')[1].split('k')[0])
        lay = int(filename.split('i')[1].split(
            'j')[1].split('k')[1].split('_')[0])
        dfMetaHss['I'].iloc[i] = row
        dfMetaHss['J'].iloc[i] = col
        dfMetaHss['K'].iloc[i] = lay

    # Get unique layer in HSS files
    layers_in_hss = dfMetaHss.K.unique()

    # total_hss_mass = 0
    df_final = pd.DataFrame()
    df2 = pd.DataFrame(columns=['Name', 'Mass'])  # pci/d or ug/d
    df2['Name'] = dfMetaHss.Name
    for i, file in enumerate(dfMetaHss.Name):
        # col name
        col_name = file[:-8]
        # print(file)

        df0 = pd.read_csv(file, delim_whitespace=True, names=[
            "time", "unknown", col_name])  # 3 columns
        # What is the unit, is it pCi?
        # print(f'size df0 = {df0.shape} /n')
        # df0 = df0[df0['time'].diff() > 0.01]
        df0['time_offset'] = df0['time'].diff()
        df0['mass_shift'] = df0[col_name].shift(periods=1)
        df0['mass_ave'] = (df0[col_name] + df0['mass_shift'])/2
        # df0['Mass_rate_shift'] = df0[col_name].shift(periods=1)
        # df0['Mass'] = df0[col_name]*df0['time_offset']
        df0['Mass'] = df0['mass_ave']*df0['time_offset']
        # if i == 0:  # Checking a cell/file
        #    print(f'df0: \n {df0}')
        #    print(f'Total mass at cell {file} = {df0["Mass"].sum()} \n')
        #
        df2['Mass'].iloc[i] = df0['Mass'].sum()/unit_conv_coef

    return df2


def post_process_df(work_dir, list_vars, dfMass, sce, model):
    today_dt = dt.datetime.today().strftime('%Y-%m-%d')
    out2 = []
    for j, var in enumerate(list_vars):
        mass = dfMass['Mass'].loc[dfMass['COPC'] == var].sum()
        out2.append([model, var, mass])
        # print(out2)

    # Get summary table
    dfMassSum = pd.DataFrame(
        data=out2, columns=['Model', 'COPC', 'Mass'])
    ofile = f'{work_dir}/output/summary_total_mass_PL_files_{sce}_{model}_{today_dt}.csv'
    dfMassSum.to_csv(
        ofile, index=False)
    print(f'Saved {ofile}\n')


def func_process_hss_dat_files(sce, list_vars, path_to_hss_dir, work_dir):
    today_dt = dt.datetime.today().strftime('%Y-%m-%d')
    # [1] for CA Rev1 Base case ---------------------------------------------------
    #list_sce = ['HSSMCAREV1']  #
    # ver = 'v1.0'
    # List of COCs for CA Rev.1 (for hss dat files)
    # list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
    #             'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']

    # [2] for  CA Rev1 Recharge Sensitivity Case ----------------------------------
    #list_sce = ['HSSMCAREV1RCH']  #
    # list_vars = ['i129', 'trit', 'tc99']

    # list_years = [1, 3, 5, 10, 15, 30, 50, 75]

    df_final = pd.DataFrame(columns=['COPC', 'Total Mass'])
    df_final['COPC'] = list_vars
    for i, var in enumerate(list_vars):
        # Get unit and conv. coefficient for var
        # unit, unit_conv_coef = var_unit(var)
        # unit, unit_conv_coef = 'pCi', 1e12  # CA, radionuclides

        # [1] read hss_files ------------------------------------------------
        path_to_hss_dir2 = os.path.join(path_to_hss_dir, var)

        # [2] to CA Rev base case
        # path_to_hss_dir = os.path.join('/',
        #                               'workspace2', 'hpham', '041_CA_rev1', 'scripts', 'input', sce, 'v1.0', 'data', var)
        # path_to_hss_dir = os.path.join('/',
        #                               'workspace2', 'hpham', '041_CA_rev1', 'scripts', 'input', sce, 'v1.0', 'data', var)

        # [3] to recharge sensitivy case of the CA Rev1 base case.
        # path_to_hss_dir = os.path.join('/',
        #                               'workspace2', 'hpham',  '041_CA_rev1', 'CA_R1_Rec_Sen', 'zip_files',  sce, 'v1.0', 'data', var)

        # CIEACM1HSSMPF (hist, ok) = 78.89608357518817 Ci
        # CIEACM1HSSM (pred)       = 371.80498784914033 Ci
        # CIEMAHSSMPF (hist)       = 78.9531122959637 Ci
        # CIEMAHSSM (pred)         = 314.586415736945 Ci
        print(f'path_to_hss_dir={path_to_hss_dir2}\n')
        dfHssTs = read_hss_files(path_to_hss_dir2, var)

        #
        df_final['Total Mass'].iloc[i] = dfHssTs["Mass"].sum()
        # print(f'dfHssTs = \n{dfHssTs} \n')
        # print(
        #    f'Total mass of {var} = {dfHssTs["Mass"].sum()} Ci\n')
        # print(f'Total mass = {df_hss_total_at_time["Mass"].sum()} Ci\n')
        # dfHssTs.to_csv(
        #    f'{work_dir}/output/total_mass_{sce}_{var}.csv')

    #
    ofile = f'{work_dir}/output/summary_total_mass_hssm_datfiles_{sce}_ver_{today_dt}.csv'
    df_final.to_csv(ofile)
    print(df_final)
    print(f'Saved {ofile}')


def func_process_VZ_facet():
    # =====================================================================
    # [2] Read VZ Facet files ---------------------------------------------
    # =====================================================================
    run_sce = 'VZ2SRIv1.1'  # VZSRIREV1 or VZ2SRIv1.1
    if run_sce == 'VZSRIREV1':
        list_sce = ['VZSRIREV1']  # CA
        ver = 'v1.1'
        # List of COCs for CA Rev.1
        list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
                     'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']
    elif run_sce == 'VZ2SRIv1.1':
        list_sce = ['VZ2SRIv1.1']  # CA
        ver = 'v1.1'
        # List of COCs for CA Rev.1
        list_vars = ['c-14', 'h-3', 'i-129']

        # df_final['No.'] = range(1, 1000, 1)  # Create 100 empty rows
    for sce in list_sce:
        out = []
        #
        path_to_data = f'/workspace2/hpham/041_CA_rev1/scripts/input/{sce}/v1.1/data/*/'
        list_dir = glob(path_to_data)
        for i, path2data in enumerate(list_dir):
            print(path2data.split())
            site = path2data.split('/')[9]
            print(site)
            # Get list of csv files in this directory path2data
            # os.chdir(path2data)
            # list_ifiles = glob(path2data+"*.csv")
            for j, var in enumerate(list_vars):
                # var = csv_file.split('/')[-1].split('-')[1]
                # print(var)

                csv_file = f'{path2data}/{site}-{var}-bot.csv'
                # print(csv_file)
                if os.path.isfile(csv_file):
                    df = pd.read_csv(csv_file, skiprows=10, header=None)
                    idx = list(range(2, df.shape[1], 2))  # [0] +
                    df2 = df.iloc[:, idx]
                    mass = df2.iloc[-1].sum()
                    #
                    out.append([sce, site, var, mass])
                    # Total mass, for example atrenches-tc-99-bot.csv = 0.05544102878386096 Ci

        #
        dfMass = pd.DataFrame(
            data=out, columns=['Sce', 'Site', 'COPC', 'Mass'])
        # dfMass.to_csv(
        #    f'{work_dir}/output/total_mass_VZ_facet_{sce}.csv')

        # post process --------------------------------------------------------
        # dfMass = pd.read_csv(
        #    f'{work_dir}/output/total_mass_VZ_facet_{sce}.csv')
        out2 = []
        for j, var in enumerate(list_vars):
            mass = dfMass['Mass'].loc[dfMass['COPC'] == var].sum()
            out2.append([sce, var, mass])
            # print(out2)

        # Get summary table
        dfMassSum = pd.DataFrame(
            data=out2, columns=['Sce', 'COPC', 'Mass'])
        dfMassSum.to_csv(
            f'/workspace2/hpham/041_CA_rev1/scripts/output/summary_total_mass_VZ_facet_{sce}.csv')


def func_process_PL_files(sce, list_vars, path2data, work_dir):
    # =====================================================================
    # [3] Read Pask Leak (PL) files ---------------------------------------
    # =====================================================================
    #work_dir = f'/workspace2/hpham/041_CA_rev1/scripts/'
    #sce = 'PAPL2SZ'
    #ver = 'v1.1'
    cutoff_year = 2018
    # cutoff_year2 = 3070 # CIE
    cutoff_year2 = 12070  # CA 12070
    print(f'Warning: Using cutoff years 2018 and 12070!!!\n')
    # list_vars = ['cn', 'cr', 'h-3', 'i-129', 'no3', 'sr-90',
    #             'tc-99', 'u-tot']
    # list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
    #             'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']

    # df_final['No.'] = range(1, 1000, 1)  # Create 100 empty rows

    out_hist, out_pred = [], []  # sce, var, mass
    unit, unit_conv_coef = 'pCi', 1e12  # CA, radionuclides
    print(f'Warning: Used unit_conv_coef={unit_conv_coef}\n')

    for i, var in enumerate(list_vars):
        #path2data = f'/workspace2/hpham/041_CA_rev1/scripts/input/{sce}/{ver}/data/{var}/'
        # /workspace2/hpham/ICF/PAPL2SZ/v1.1/data/tc-99
        # list_ifiles = glob(
        #    path2data+"WMA_C_PL*.csv") # CIE
        path2data2 = os.path.join(path2data, var)
        list_ifiles = glob(
            path2data2+"/*.csv")  # CA

        # unit, unit_conv_coef = var_unit(var) # CIE

        print(f'{var}, list_ifiles={list_ifiles}\n')
        for j, csv_file in enumerate(list_ifiles):
            # print(i, j, csv_file)
            # csv_file = f'/workspace2/hpham/ICF/PAPL2SZ/v1.1/data/cr/WMA_C_PL_cr_rates_total_area_yearly_steps.csv'
            df = pd.read_csv(csv_file, skiprows=3, header=None)
            # Historical data
            df_hist = df.loc[df.iloc[:, 0] < cutoff_year]
            # df_hist=df_hist.reset_index()

            df_pred = df.loc[df.iloc[:, 0] >= cutoff_year]
            df_pred = df_pred.loc[df_pred.iloc[:, 0] < cutoff_year2]
            df_pred = df_pred.reset_index(drop=True)

            # if i == 1:
            #    print(df_hist.shape, df_hist.head())
            #    print(df_pred.shape, df_pred.head())

            # Hist dataframe
            # get data from 2nd column
            idx = list(range(1, df.shape[1], 1))
            df_hist2 = df_hist.iloc[:, idx]
            idx2 = list(range(0, df_hist2.shape[0], 2))  # [0] +
            df_hist3 = df_hist2.iloc[idx2, :]
            mass_hist = df_hist3.sum().sum()/unit_conv_coef
            out_hist.append([sce, csv_file, var, mass_hist])

            # Predictive dataframe
            df_pred2 = df_pred.iloc[:, idx]
            idx2 = list(range(0, df_pred2.shape[0], 2))  # [0] +
            df_pred3 = df_pred2.iloc[idx2, :]
            mass_pred = df_pred3.sum().sum()/unit_conv_coef
            out_pred.append([sce, csv_file, var, mass_pred])

            #
            # Total mass, for example atrenches-tc-99-bot.csv = 0.05544102878386096 Ci

        # hist
        dfMass_hist = pd.DataFrame(
            data=out_hist, columns=['Sce', 'Site', 'COPC', 'Mass'])
        # dfMass_hist.to_csv(
        #    f'{work_dir}/output/total_mass_PL_files_hist.csv')
        # pred
        dfMass_pred = pd.DataFrame(
            data=out_pred, columns=['Sce', 'Site', 'COPC', 'Mass'])
        # dfMass_pred.to_csv(
        #    f'{work_dir}/output/total_mass_PL_files_pred.csv')

        # postprocessing each scenario ------------------------------------
        post_process_df(work_dir, list_vars, dfMass_hist, sce, 'hist')
        post_process_df(work_dir, list_vars, dfMass_pred, sce, 'pred')


def func_process_sri_files(list_vars, sce, path2data, work_dir):
    # =====================================================================
    # [4] Read SRI files --------------------------------------------------
    # =====================================================================
    # opt_task = 'VZ2SRIv1.1'
    today_dt = dt.datetime.today().strftime('%Y-%m-%d')
    '''
    if opt_task == 'CIE_MA':
        list_sce = ['CIEMASRI2018', 'CIEMASRI3070']  #
        ver = 'v1.0'
        # list_vars = ['cn', 'cr', 'h-3', 'i-129', 'no3', 'sr-90',
        #             'tc-99', 'u']  # CIE MA
        # List of COCs for CA Rev.1
        list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
                     'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']

    elif opt_task == 'CA_Rev1':
        list_sce = ['SRICAREV1']  #
        ver = 'v1.0'
        # List of COCs for CA Rev.1
        list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
                     'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']

    elif opt_task == 'SRICAREV1MOD':
        list_sce = ['SRICAREV1MOD']  #
        ver = 'v1.0'
        # List of COCs for CA Rev.1
        list_vars = ['i-129', 'h-3', 'tc-99']

    elif opt_task == 'CA_Rev1_RCH_SEN':
        list_sce = ['SRICAREV1RCH']  #
        ver = 'v1.0'
        # List of COCs for CA Rev.1
        list_vars = ['i-129', 'h-3', 'tc-99']
    '''
    # df_final['No.'] = range(1, 1000, 1)  # Create 100 empty rows
    # for sce in list_sce:
    out = []  # sce, var, mass
    for i, var in enumerate(list_vars):
        '''
        if (opt_task == 'CIE_MA' or opt_task == 'CA_Rev1'):
            path2data = f'/workspace2/hpham/041_CA_rev1/scripts/input/{sce}/{ver}/data/{var}/'
        elif (opt_task == 'CA_Rev1_RCH_SEN' or opt_task == 'SRICAREV1MOD'):
            path2data = f'/workspace2/hpham/041_CA_rev1/CA_R1_Rec_Sen/zip_files/{sce}/{ver}/data/{var}/'
            print(f'path2data={path2data}\n')
        '''
        path2data2 = os.path.join(path2data, var)
        print(f'path2data2={path2data2}\n')

        list_ifiles = glob(path2data2+"/*.csv")
        print(f'list_ifiles={list_ifiles}\n')

        # unit, unit_conv_coef = var_unit(var) # CIE
        unit, unit_conv_coef = 'pCi', 1e12  # CA, radionuclides
        # csv_file=f'/workspace2/hpham/ICF/CIEMASRI3070/v1.0/data/tc-99/afarms-tc-99-bot_yearly_steps.csv'
        print(f'Warning: Using unit_conv_coef={unit_conv_coef}\n')
        for j, csv_file in enumerate(list_ifiles):
            df = pd.read_csv(csv_file, skiprows=3, header=None)
            idx = list(range(1, df.shape[1], 1))  # [0] +
            df2 = df.iloc[:, idx]
            # Get every other cell
            idx2 = list(range(0, df2.shape[0], 2))  # [0] +
            df3 = df2.iloc[idx2, :]
            mass = df3.sum().sum()/unit_conv_coef
            #
            out.append([sce, csv_file, var, mass])
            # Total mass, for example atrenches-tc-99-bot.csv = 0.05544102878386096 Ci

        #
        dfMass = pd.DataFrame(
            data=out, columns=['Sce', 'Site', 'COPC', 'Mass'])
        # dfMass.to_csv(
        #    f'{work_dir}/output/total_mass_sri_file_{sce}.csv')

    # postprocessing each scenario ------------------------------------
    # dfMass = pd.read_csv(
    #    f'{work_dir}/output/total_mass_sri_file_{sce}.csv')
    out2 = []
    for j, var in enumerate(list_vars):
        mass = dfMass['Mass'].loc[dfMass['COPC'] == var].sum()
        out2.append([sce, var, mass])
        # print(out2)

    # Get summary table
    dfMassSum = pd.DataFrame(
        data=out2, columns=['Sce', 'COPC', 'Mass'])
    ofile = f'{work_dir}/output/summary_total_mass_sri_file_{sce}_{today_dt}.csv'
    dfMassSum.to_csv(ofile)
    print(dfMassSum)
    print(f'Save {ofile}\n')


def func_map_mass_distribution(run_sce):
    today_dt = dt.datetime.today().strftime('%Y-%m-%d')
    print('Running opt_map_hss_mass_distribution')
    dic_levels = {'c-14': [0, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1, 5e12],
                  'cl36': [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e12],
                  'cl-36': [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e12],
                  'trit': [0, 5e-4, 1e-3, 5e-3, 1e-2, 2e-2, 1e12],
                  'h-3': [0, 5e-4, 1e-3, 5e-3, 1e-2, 2e-2, 1e12],
                  'i129': [0, 1e-7, 1e-6, 2e-6, 5e-6, 1e-5, 1e12],
                  'i-129': [0, 1e-7, 1e-6, 2e-6, 5e-6, 1e-5, 1e12],
                  'n237': [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e12],
                  'np-237': [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e12],
                  'r226': [0, 1e-10, 5e-10, 1e-9, 2e-9, 3e-9, 1e12],
                  'ra-226': [0, 1e-10, 5e-10, 1e-9, 2e-9, 3e-9, 1e12],
                  'sr90': [0, 1e-8, 5e-8, 1e-7, 2e-7, 4e-7, 1e12],
                  'sr-90': [0, 1e-8, 5e-8, 1e-7, 2e-7, 4e-7, 1e12],
                  'tc99': [0, 1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 1e12],
                  'tc-99': [0, 1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 1e12],
                  't230': [0, 1e-11, 3e-11, 9e-11, 1e-10, 4e-10, 1e12],
                  'th-230': [0, 1e-11, 3e-11, 9e-11, 1e-10, 4e-10, 1e12],
                  'u232': [0, 1e-11, 5e-11, 1e-10, 5e-10, 1e-9, 1e12],
                  'u-232': [0, 1e-11, 5e-11, 1e-10, 5e-10, 1e-9, 1e12],
                  'u233': [0, 1e-7, 2e-7, 5e-7, 1e-6, 5e-6, 1e12],
                  'u-233': [0, 1e-7, 2e-7, 5e-7, 1e-6, 5e-6, 1e12],
                  'u234': [0, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e12],
                  'u-234': [0, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e12],
                  'u235': [0, 1e-8, 5e-8, 1e-7, 1e-6, 2e-6, 1e12],
                  'u-235': [0, 1e-8, 5e-8, 1e-7, 1e-6, 2e-6, 1e12],
                  'u236': [0, 1e-8, 1e-7, 2e-7, 5e-7, 1e-6, 1e12],
                  'u-236': [0, 1e-8, 1e-7, 2e-7, 5e-7, 1e-6, 1e12],
                  'u238': [0, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 1e12],
                  'u-238': [0, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 1e12]}
    list_vars = ['c-14', 'cl36', 'trit', 'i129', 'n237', 'r226',  'sr90',  'tc99',  't230',
                 'u232', 'u233', 'u234', 'u235', 'u236', 'u238']

    if run_sce == 'CA_Rev0':
        # CA Rev0 files: --------------------------------------------------
        # This file is from ICF/Prod
        # path_to_hss = os.path.join('/',
        #                           'workspace2', 'hpham', '041_CA_rev1', 'scripts', 'input', 'HSSMCA4CARev0', 'data')

        # This is actual hssm files used in CA Rev0, in nsbs.tar
        path_to_hss = os.path.join('/',
                                   'workspace2', 'hpham', '041_CA_rev1', 'scripts', 'input', 'HSSMRev0_TB', 'hssm', 'nonStep')

        # list_vars = ['c-14', 'cl36', 'trit', 'i129', 'n237', 'r226',  'sr90',  'tc99',  't230',
        #             'u232', 'u233', 'u234', 'u235', 'u236', 'u238']

    elif run_sce == 'CA_Rev1':
        # CA Rev1 files: ------------------------------------------------------
        path_to_hss = os.path.join('/',
                                   'workspace2', 'hpham', '041_CA_rev1', 'scripts',
                                   'input', 'HSSMCAREV1', 'v1.0', 'data_rev')
        # list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
        #             'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238']
        list_vars = ['c-14', 'cl36', 'trit', 'i129', 'n237', 'r226',  'sr90',  'tc99',  't230',
                     'u232', 'u233', 'u234', 'u235', 'u236', 'u238']

        # to compare with REC SENS run ----------------------------------------
        # list_vars = ['i129', 'trit', 'tc99']
        work_dir2 = os.path.join('/', 'workspace2', 'hpham',
                                 '041_CA_rev1')
        # end comparison ------------------------------------------------------
    elif run_sce == 'CA_Rev0_GAIA':
        # CA Rev1 files: ------------------------------------------------------
        path_to_hss = os.path.join('/',
                                   'workspace2', 'hpham', '041_CA_rev1',
                                   'zip_files', 'CArev0_hssm_20211117', 'hssm')
        # list_vars = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
        #             'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238']
    elif run_sce == 'CA_r1_REC_SEN':
        path_to_hss = os.path.join('/',
                                   'workspace2', 'hpham', '041_CA_rev1',
                                   'CA_R1_Rec_Sen', 'zip_files', 'HSSMCAREV1RCH', 'v1.0', 'data')
        list_vars = ['i129', 'trit', 'tc99']
        # dir to save ouput such as csv, txt, png files
        work_dir2 = os.path.join('/', 'workspace2', 'hpham',
                                 '041_CA_rev1', 'CA_R1_Rec_Sen')
    # /workspace2/hpham/041_CA_rev1/CA_R1_Rec_Sen/zip_files/HSSMCAREV1RCH/v1.0/data
    # 'ECF-HANFORD-21-0043_R0' and 'ECF-HANFORD-21-0044_R0'
    # ecf = 'ECF-HANFORD-21-0043_R0' # CIE
    # ecf = 'CA Rev1'
    # list_vars = ['cn', 'cr', 'trit', 'i129', 'no3', 'sr90',
    #             'tc99', 'u']  # 'ECF-HANFORD-21-0043_R0'
    # list_vars = ['i129', 'tc99']  # 'ECF-HANFORD-21-0043_R0'

    # Read P2RGWM Version 8.3 Grid Information
    ifile_grid = os.path.join(work_dir, 'input', 'grid_274_geo_coor.csv')
    dfGrid = pd.read_csv(ifile_grid)

    # Read hss dat files --------------------------------------------------
    # sce, ver = 'CIEMAHSSM', 'v1.0'  # CIE MA
    # sce, ver = 'CIENFAHSSM', 'v1.2'  # CIE NFA, for QA this script
    # sce, ver = 'CIEACM1HSSM', 'v1.0'  # CIE ACM1, for QA this script
    # sce, ver = 'HSSMCAREV1', 'v1.0'
    for var in list_vars:
        # Add TB's HSS val to dfGrid for checking
        '''
        ifile_hss_TB = f'/workspace2/hpham/041_CA_rev1/fr_TB/mass_ref_files/{var}.ref'
        arr_TB = read_ref(ifile_hss_TB, nr, nc)
        dfGrid['HSS_TB'] = np.reshape(arr_TB, nr*nc)
        print(f'path_to_hss={path_to_hss}\n')
        '''
        #
        path_to_hss_dir = os.path.join(path_to_hss, var)
        # /workspace2/hpham/ICF/CIENFAHSSM/v1.2/data/i129
        # 041_CA_rev1/scripts/input/HSSMCAREV1/v1.0

        df = read_hss_files(path_to_hss_dir, var)
        ROW, COL, LAY = [], [], []
        for i in range(df.shape[0]):
            ROW.append(
                int(re.split('i|j|k|_', df.Name.iloc[i])[1]))  # irow
            COL.append(
                int(re.split('i|j|k|_', df.Name.iloc[i])[2]))  # irow
            LAY.append(
                int(re.split('i|j|k|_', df.Name.iloc[i])[3]))  # irow
        #
        df['ROW'], df['COL'], df['LAY'] = ROW, COL, LAY
        # Save df for checking
        # df.to_csv(os.path.join(work_dir, 'output', 'check_df.csv'))

        # Add cell area to df -------------------------------------------------
        df_final = pd.merge(df, dfGrid, how="left", on=['ROW', 'COL'])
        df_final['Mass_per_m2'] = df_final['Mass']/df_final['area']

        # Create array
        arr = np.empty([nr, nc], dtype='float')

        for k in range(df_final.shape[0]):
            i = df_final['I'].iloc[k] - 1
            j = df_final['J'].iloc[k] - 1
            arr[i, j] = df_final['Mass_per_m2'].iloc[k]
        #
        df_final = df_final[df_final['Mass_per_m2'] != 0]
        ofile_csv2 = os.path.join(
            work_dir, 'output', 'csv', 'total_activity', f'check_total_activity_{run_sce}_{var}_{today_dt}.csv')
        df_final.to_csv(ofile_csv2)
        print(f'Saved {ofile_csv2}\n')

        # [1] Save to shp file format -----------------------------------------
        ofile = os.path.join(work_dir2, 'scripts',
                             'output', 'shp', f'{run_sce}_avg_{var}.shp')
        # arr2shp(work_dir2, arr, ofile)
        # print(f'Saved {ofile}\n')

        # [2] save to png file ------------------------------------------------
        # out_dir = os.path.join('/', 'workspace2', 'hpham',
        #                       '041_CA_rev1', 'scripts')  # CIE MA
        ofile_png = f'{work_dir2}/scripts/output/png/hss_total_mass/hss_total_mass_{var}_{run_sce}_{today_dt}.png'
        ptitle = f'Total estimated activity arriving at water table for {var}, {run_sce}'
        colors = ['#63b8ff', '#bfefff', '#bab36a',
                  '#ffff00', '#f8b186', '#f70fe8']  # hss

        map_type = 'hss'
        # levels = [1e2, 1e3, 1e4, 1e5, 1e6, 1e9, 1e10, 1e11]
        arr = np.ma.masked_less_equal(arr, 1e-20)
        vmin, vmax = np.nanmin(arr), np.nanmax(arr)
        print(vmin, vmax)
        generate_map1(arr, ofile_png,  # DIFF_percent
                      f'{work_dir2}/scripts/', ptitle, dic_levels[var], colors, map_type, ca_cie=True)
        print(f'Saved {ofile_png}\n')

        # [3] Save arr to text file -------------------------------------------
        ofile_arr = os.path.join(work_dir2, 'scripts',
                                 'output', 'txt', f'{run_sce}_avg_{var}_{today_dt}.txt')
        np.savetxt(ofile_arr, arr, fmt='%.15f', delimiter=',')
        print(f'Saved {ofile_arr}\n')


if __name__ == "__main__":
    # [1] Load input file =========================================================
    #ifile = f'input/input_parameters.csv'
    # ifile = f'input/input_CA_R1_INV_ST.csv'
    ifile = f'input/input_CA_R1_NSMC_UCN.csv'

    #
    dfin = pd.read_csv(ifile)
    dfin = dfin.set_index('var')
    sce = dfin['name'].loc['sce']
    list_vars = re.split(', ', dfin['name'].loc['list_vars'])
    work_dir = dfin['name'].loc['pth2scpt']  # full path to this cript
    # full path to this cript
    nlay, nr, nc = re.split(', ', dfin['name'].loc['nlay_nr_ncol'])
    nlay, nr, nc = int(nlay), int(nr), int(nc)

    #
    today_dt = dt.datetime.today().strftime('%Y-%m-%d')
    # work_dir = os.path.join('/', 'workspace2', 'hpham',
    #                        '034_CIE_MA', 'scripts')
    # work_dir = os.path.join('/', 'workspace2', 'hpham',
    #                        '041_CA_rev1', 'scripts')

    # nlay, nr, nc = 7, 201, 274
    # list_sce = ['CIEMAHSSMPF', 'CIEMAHSSM']  # MA scenarios
    # list_sce = ['CIEACM1HSSMPF', 'CIEACM1HSSM']  # ACM1 NFA scenarios 'CIEACM1HSSMPF', 'CIEACM1HSSM'
    # list_sce = ['CIEACM1MAHSSM']  # CIE ACM1 MA
    # list_vars = ['cn', 'cr', 'trit', 'i129', 'no3', 'sr90',
    #             'tc99', 'u']
    # list_vars = ['i129', 'tc99']

    # [0] Choose run option ---------------------------------------------------
    opt_hss_dat_files = eval(dfin['name'].loc['opt_hss_dat_files'])
    opt_VZ_facet_files = eval(dfin['name'].loc['opt_VZ_facet_files'])
    opt_PL_files = eval(dfin['name'].loc['opt_PL_files'])
    opt_sri_files = eval(dfin['name'].loc['opt_sri_files'])
    opt_P2RGWM_mas_files = eval(dfin['name'].loc['opt_P2RGWM_mas_files'])
    # generate png, shp, txt, csv of mass/activity distribution
    opt_map_hss_mass_distribution = eval(
        dfin['name'].loc['opt_map_hss_mass_distribution'])
    # png shp map hss diff rev0 vs. rev1
    opt_check_hss_diff = eval(dfin['name'].loc['opt_check_hss_diff'])
    # combine png rev0 and rev1 ...
    opt_comb_hss_mas_files = eval(dfin['name'].loc['opt_comb_hss_mas_files'])
    # Compare two ucn files
    opt_compare_2ucn = eval(dfin['name'].loc['opt_compare_2ucn'])

    # ==========================================================================
    # END OF READING INPUT FILE ================================================
    # ==========================================================================

    # Create folder to save output
    create_new_dir('output')
    # [1] Read and process hss dat files --------------------------------------
    if opt_hss_dat_files:
        path_to_hss_dir = dfin['path2files'].loc['opt_hss_dat_files']
        func_process_hss_dat_files(sce, list_vars, path_to_hss_dir, work_dir)
        os.chdir(work_dir)

    # [2] Read VZ Facet files -------------------------------------------------
    if opt_VZ_facet_files:
        func_process_VZ_facet(sce, list_vars, path_to_hss_dir, work_dir)

    # [3] Read Pask Leak (PL) files -------------------------------------------
    if opt_PL_files:
        path_to_PL_files = dfin['path2files'].loc['opt_PL_files']
        func_process_PL_files(sce, list_vars, path_to_PL_files, work_dir)

    # [4] Read SRI files ------------------------------------------------------
    if opt_sri_files:
        # choose among 'CIE_MA' or 'CA_Rev1' or 'CA_Rev1_RCH_SEN' or 'VZ2SRIv1.1'
        # SRICAREV1MOD: HP removed some files that are already in CA_Rev1_RCH_SEN
        #opt_task = 'SRICAREV1MOD'
        path2data = dfin['path2files'].loc['opt_sri_files']
        print(f'\n\npath2data={path2data}\n')
        func_process_sri_files(list_vars, sce, path2data, work_dir)

    # [5] Read P2RGWM.mas files -----------------------------------------------
    if opt_P2RGWM_mas_files:
        # 'acm1': ACM1 scenario, 'base': MA scenarios in this;
        '''
        ecf = 'CA_r1_base'  # 'CA_r1_base' or 'CA_r1_REC_SEN' or 'CA_r1_435_SEN'
        # work_dir = os.path.join('/', 'workspace2', 'hpham',
        #                        '034_CIE_MA') # CIE MA
        # work_dir = os.path.join('/', 'workspace2', 'hpham',
        #                        '041_CA_rev1', 'rev1')  # CIE MA rev2.0a
        if ecf == 'CA_r1_base':
            work_dir = os.path.join('/', 'workspace2', 'hpham',  # CIE MA rev2.0
                                    '041_CA_rev1', 'base', 'tran', 'rev2.0_ST',  'avg')
            # dir to save output
            script_dir = os.path.join('/', 'workspace2', 'hpham',
                                      '041_CA_rev1')  #

        elif ecf == 'CA_r1_REC_SEN':
            work_dir = os.path.join('/', 'workspace2', 'hpham',  # CIE MA rev2.0
                                    '041_CA_rev1', 'CA_R1_Rec_Sen', 'model_runs', 'tran', 'nsbs', 'avg')
            list_vars = ['i129', 'trit', 'tc99']
            script_dir = os.path.join('/', 'workspace2', 'hpham',
                                      '041_CA_rev1', 'CA_R1_Rec_Sen')  #
        elif ecf == 'CA_r1_435_SEN':
            work_dir = os.path.join('/', 'workspace2', 'hpham',  # CIE MA rev2.0
                                    '041_CA_rev1', 'CA_R1_435_Sen', 'model_runs', 'tran', 'nsbs', 'avg')
            list_vars = ['i129', 'tc99', 'u238']
            script_dir = os.path.join('/', 'workspace2', 'hpham',
                                      '041_CA_rev1', 'CA_R1_435_Sen')  #
        '''

        # Save table to csv
        path2masfile = dfin['path2files'].loc['opt_P2RGWM_mas_files']
        dfP2R_Mas = read_mas_files(path2masfile, list_vars)
        dfP2R_Mas['Mass_Ci'] = dfP2R_Mas['Mass']/1e12
        print(dfP2R_Mas)
        ofile = f'{work_dir}/output/check_dfP2R_Mas_Table_{sce}_HP_{today_dt}.csv'
        dfP2R_Mas.to_csv(ofile, index=False)
        print(f'Saved {ofile}\n')

    # [6] Map spatial distribution of total mass/activity arriving at water table
    if opt_map_hss_mass_distribution:
        '''
        This generates:
           check_total_activity_{run_sce}_{var}.csv # For checking
           {run_sce}_avg_{var}.shp # for mapping
           {run_sce}_avg_{var}.txt # array to check mass diff
        '''
        # CA_Rev0 or CA_Rev1 or CA_Rev0_GAIA or CA_r1_REC_SEN
        run_sce = 'CA_Rev1'  # CA_Rev1 (means base) or CA_R1_Sen
        func_map_mass_distribution(run_sce)
    #
    # [7] Plot hss differences
    # run firstly opt_map_hss_mass_distribution
    if opt_check_hss_diff:
        list_vars2 = ['c-14', 'cl-36', 'h-3', 'i-129', 'np-237', 'ra-226',  'sr-90',  'tc-99',  'th-230',
                      'u-232', 'u-233', 'u-234', 'u-235', 'u-236', 'u-238', 're-187']
        for k, coc in enumerate(list_vars):
            # for k, coc in enumerate(['c-14']):
            ifile0 = f'output/txt/CA_Rev0_avg_{coc}.txt'
            v0 = np.loadtxt(ifile0, delimiter=',')

            coc1 = list_vars2[k]
            ifile1 = f'output/txt/CA_Rev1_avg_{coc1}.txt'
            v1 = np.loadtxt(ifile1, delimiter=',')

            #
            v0 = np.multiply(v0, 1e12)  # convert to PCi
            v1 = np.multiply(v1, 1e12)  # convert to PCi
            DIFF = v1-v0  # PCi

            # DIFF = v1-v0
            DIFF = np.ma.masked_equal(DIFF, 0.0)  # Mask 0 cells
            min_ = np.nanmin(DIFF)
            max_ = np.nanmax(DIFF)
            mean_ = np.nanmean(DIFF)
            print(min_, max_, mean_)

            #
            nr, nc = DIFF.shape
            print(nr, nc)
            val0 = np.reshape(v0, nr*nc)
            val1 = np.reshape(v1, nr*nc)
            val = np.reshape(DIFF, nr*nc)

            # use model grid shapefile geometry
            # model grid shapefile
            gridShp = os.path.join('/workspace2/hpham/041_CA_rev1', 'shp_files',
                                   'grid_274_geo_rc.shp')

            gdf = gpd.read_file(gridShp)
            df = pd.DataFrame()
            df['row'] = gdf['row']
            df['column'] = gdf['column']
            df['rev0_pCi_per_m2'] = val0
            df['rev1_pCi_per_m2'] = val1
            #
            df['Diff'] = val
            df['Percent_Diff'] = 100 * \
                (df['rev1_pCi_per_m2']-df['rev0_pCi_per_m2']) / \
                df['rev0_pCi_per_m2']

            #
            work_dir = os.path.join('/', 'workspace2', 'hpham',
                                    '041_CA_rev1', 'scripts')  # CIE MA

            # [2] to png file format ----------------------------------------------
            ofile_png = f'{work_dir}/output/png/check_hss_diff_{coc}.png'
            # fig, ax = plt.subplots()
            ptitle = f'Difference in total activity (pCi/m2) for {coc}. DIFF=abs(Rev1-Rev0)'
            colors = ['#63b8ff', '#bfefff', '#bab36a',
                      '#ffff00', '#f8b186', '#f70fe8']  # hss

            map_type = 'hss'
            # generate_map_CIE(ax, fig, DIFF, work_dir, cscale,
            #              ptitle, ca_cie=False)
            print(df['Percent_Diff'].shape)
            DIFF_percent = np.reshape(df['Percent_Diff'].to_numpy(), (nr, nc))
            # levels = [1e2, 1e3, 1e4, 1e5, 1e6, 1e9, 1e10, 1e11]
            # generate_map1(np.abs(DIFF), ofile_png,  # DIFF_percent
            #              work_dir, ptitle, levels, colors, map_type, ca_cie=True)
            ofile_png = f'{work_dir}/output/png/check_hss_diff_percent/check_hss_diff_percent_{coc}.png'
            levels = [0, 0.1, 0.5, 1, 2.5, 25, 50, 100]
            generate_map1(np.abs(DIFF_percent), ofile_png,  # DIFF_percent
                          work_dir, ptitle, levels, colors, map_type, ca_cie=True)

            # [3] to shp file format ----------------------------------------------
            ofile = f'{work_dir}/output/shp/check_hss_diff_percent/check_hss_diff_percent_{coc}.shp'
            arr2shp(DIFF_percent, DIFF, ofile)
            print(f'Save {ofile}\n')

            #
            # fig.savefig(ofile_png, dpi=150, transparent=False, bbox_inches='tight')
            # print(f'The figure was saved at: {ofile_png}')
            plt.clf()
            plt.close('all')
            # [4] to csv for checking
            df = df[df['Diff'] != np.nan]
            ofile_csv = f'{work_dir}/output/csv/check_hss_diff_{coc}.csv'
            df.to_csv(ofile_csv, index=False)
    # [8]
    if opt_comb_hss_mas_files:
        '''
        merge csv files
        '''
        coc = 'c-14'
        ifile0 = f'output/csv/total_activity/check_total_activity_CA_Rev0_{coc}.csv'
        ifile1 = f'output/csv/total_activity/check_total_activity_CA_Rev1_{coc}.csv'
        df0 = pd.read_csv(ifile0)
        df0.rename(columns={"Mass": "MasRev0"}, inplace=True)

        df1 = pd.read_csv(ifile1)

        df1 = df1[['Name', 'Mass']]
        df1.rename(columns={"Mass": "MasRev1"}, inplace=True)

        df_merge = df0.merge(df1, how='left', on='Name')
        df_merge = df_merge[['Name', 'MasRev0', 'MasRev1', 'ROW', 'COL', 'LAY', 'FID', 'delx',
                             'dely', 'area', 'I', 'J', 'IBND', 'X', 'Y']]
        df_merge['Diff_r1_r0'] = df_merge['MasRev1']-df_merge['MasRev0']
        df_merge['Abs_Diff_pCi'] = 1e12*df_merge['Diff_r1_r0'].abs()
        df_merge.to_csv(
            f'output/csv/total_activity/diff_{coc}.csv', index=False)

    # =========================================================================
    # [9] Ref file to dataframe -----------------------------------------------
    # =========================================================================
    '''
    Ref files of total activitity arriving at water table, calculated by TB.
    This script reshapes ref files to arrays of row, col and val
    '''
    # list_vars = ['c-14', 'cl36', 'trit', 'i129', 'n237', 'r226',  'sr90',  'tc99',  't230',
    #             'u232', 'u233', 'u234', 'u235', 'u236', 'u238']
    opt_process_TB_ref_files = False
    if opt_process_TB_ref_files:
        for k, coc in enumerate(list_vars):
            # Ref file of hss mass from TB (email 11/16/2021)
            #
            ifile_hss_TB = f'/workspace2/hpham/041_CA_rev1/fr_TB/mass_ref_files/{coc}.ref'
            arr = read_ref(ifile_hss_TB, nr, nc)

            # model grid shapefile
            gridShp = os.path.join('/workspace2/hpham/041_CA_rev1', 'shp_files',
                                   'grid_274_geo_rc.shp')

            gdf = gpd.read_file(gridShp)
            df = pd.DataFrame()
            df['row'] = gdf['row']
            df['column'] = gdf['column']
            df['val'] = np.reshape(arr, nr*nc)

            # save file
            ofile = f'output/csv/total_activity_TB_ref/total_activity_TB_{coc}.csv'
            df.to_csv(ofile, index=False)
            print(f'Saved {ofile}\n')

            # Save only hss cells
    #
    #
    #
    if opt_compare_2ucn:
        path2files = re.split(', ', dfin['path2files'].loc['opt_compare_2ucn'])
        func_compare_2ucn_files(sce, list_vars, path2files)
