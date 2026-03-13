#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare groups microstrcture
@author: localadmin
"""

import os
import sys
import matplotlib.pyplot as plt
from dmri_dmrs_toolbox.dmrs.dmrsmodel import DMRSModel

subj_list = [f'sub-{i:02d}' for i in [3,4,5,6,7,8,10,11,12,14,15]]    # list of subjects to analyse [8,10,11,12,14,15]]#

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = "/media/localadmin/DATA/data/CTD/"          # path to where the data from the cohort is
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = "/home/localadmin/Software/dMRI_dMRS_toolbox/common/"  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList_CTD.xlsx'   # name of the excel file containing the metadata of the cohort
cfg['model_list']           =  ['cylinder']
cfg['metabolites']          = ['NAA+NAAG','Glu','Ins','GPC+PCho','Cr+PCr','Tau','Gln']              # metabolites for analysis


from src.dmri_dmrs_toolbox.misc.bids_structure import *
from src.dmri_dmrs_toolbox.misc.custom_functions import *

from scipy.optimize import curve_fit
import pandas as pd
import glob
import copy
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu
from pathlib import Path

output_folder = Path(cfg['data_path'])/'results'
output_folder.mkdir(parents=True, exist_ok=True)

scan_list   = pd.read_excel(os.path.join(cfg['data_path'], cfg['scan_list_name']))

def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx - 1)*np.var(x, ddof=1) + (ny - 1)*np.var(y, ddof=1)) / dof)



######## MODEL-WISE OPERATIONS ########
dmrsmodel = DMRSModel()
for model in cfg['model_list']:
    print(f'Working with {model}...')
    params = dmrsmodel.model_specs[model].param_names
    nb_params = len(params)

    df_all_data = pd.DataFrame()

    for subj in cfg['subj_list']:
        print('Working on subject ' + subj + '...')
        # Extract data for subject
        subj_data = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        data_path = cfg['data_path']
        group = subj_data['group'].iloc[0]
        if '?' in group:
            group = 'no group yet'

        # List of acquisition sessions
        sess_list = [x for x in list(subj_data['sessNo'].unique()) if not math.isnan(x)]  # clean NaNs

        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            print('Working on session ' + str(sess) + '...')

            bids_strc_analysis = create_bids_structure(subj=subj, sess=sess, datatype='dmrs', root=data_path,
                                                       folderlevel='derivatives', workingdir=cfg['analysis_foldername'],
                                                       description=model)

            for metab in cfg['metabolites']:
                fit_params_path = Path(bids_strc_analysis.get_path())/"csvs"/f'fit_parameters_{metab}_{model}.csv'
                if fit_params_path.is_file():
                    fit_parameters = pd.read_csv(fit_params_path)
                else:
                    fit_parameters = pd.DataFrame(
                        [np.nan]*nb_params,
                        columns = params
                    )
                    group = "no data"
                fit_parameters['Metabolite'] = metab
                fit_parameters["Subject"] = subj
                fit_parameters["Group"] = group
                fit_parameters["Session"] = int(sess)

                df_all_data = pd.concat([df_all_data, fit_parameters], ignore_index=True)

    ######## FOR ALL THE DATA OF THIS MODEL ########
    print("hi")
