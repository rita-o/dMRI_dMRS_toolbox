#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to analyse dMRS data

Last changed Jan 2025
@author: Rita, Malte
"""
 
import os
import sys
import matplotlib.pyplot as plt
import json
from datetime import datetime
from pathlib import Path

plt.close('all');
os.system('clear')
os.system('cls')

########################## SCRIPT CONFIGURATION (EDIT AS APPPROPRIATE) ##########################

#### DATA PATH AND SUBJECTS ####
subj_list = [f'sub-{i:02d}' for i in [3,4,5]] #,3,4,5]]    # list of subjects to analyse

cfg                         = {}
cfg['subj_list']            = subj_list
cfg['data_path']            = "/media/localadmin/DATA/data/BDL/"          # path to where the data from the cohort is
cfg['code_path']            = "/media/localadmin/DATA/data/BDL/dMRI-dMRS-Processing-Toolbox/" # path to code folder
cfg['code_path2']           = str(Path(cfg['code_path']) / 'processing_dmrs')    # path to code subfolder
cfg['toolboxes']            = "/home/localadmin/Software/"                              # path to where some toolboxes from matlab are (including MPPCA and tMPPCA)
cfg['prep_foldername']      = 'preprocessed'    # name of the preprocessed folder (keep 'preprocessed' as default)
cfg['analysis_foldername']  = 'analysis'        # name of the analysis folder (keep 'analysis' as default)
cfg['common_folder']        = str(Path(cfg['code_path'])/'common')  # path to the common folder with files needed throught the pipeline
cfg['scan_list_name']       = 'ScanList_BDL.xlsx'   # name of the excel file containing the metadata of the cohort

#### ADD CODE PATH ####
sys.path.append(cfg['code_path'])
sys.path.append(cfg['code_path2'])

import importlib
from bids_structure import *
from custom_functions import *

#### DMRS PREPROCESSING CONFIG ####

cfg['LC_model']             = os.path.join(cfg['toolboxes'], 'LCModel','binaries','linux')  # path to LC model executable
cfg['MATLAB_Runtime']       = "/home/localadmin/Software/Matlab_runtime/R2025a" # path to MATLAB_Runtime. It needs to be R2025a: https://ch.mathworks.com/products/compiler/matlab-runtime.html
cfg['basis_set']            = os.path.join(cfg['common_folder'], 'mrs_basis_sets','Basis_Set_dSPECIAL_differentTM')     # path to where the basis set are
cfg['models']               = ["dti","stick", "dki","cylinder","cylinder_sphere","stick_sphere"]    # models used for fitting
cfg['metabolites']          = ['NAA+NAAG','Glu','Ins','GPC+PCho','Cr+PCr','Tau','Gln']              # metabolites for analysis
# for subject 2, tCho is not working, use this reduced list of metabs instead
#cfg['metabolites']          = ['NAA+NAAG','Glu','Ins','Cr+PCr','Tau','Gln']              # metabolites for analysis
cfg['redo_processing']      = 1  # 1 to remove previous file and redo all processing; 0 to process only missing steps

#### SAVE CONFIG FILE ####
with open(Path(cfg['data_path']) / '.config_mrs.json', 'w') as f:
    json.dump(cfg, f)

#### STEP 1. Process and quantify bruker data
from Step1_preproc import *
#Step1_preproc(cfg)

#### STEP 2. Fitting of data (needs SwissKnife environment)
from Step2_fitting import *

analysis_dir = Path(cfg['data_path'])/'derivatives'/cfg['analysis_foldername']
analysis_dir.mkdir(parents=True, exist_ok=True)

Step2_fitting(cfg)

env = os.environ.copy()
env["QT_QPA_PLATFORM"] = "offscreen"
env["XDG_RUNTIME_DIR"] = "/tmp"
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
path_log = os.path.join(cfg['data_path'],'derivatives',cfg['analysis_foldername'],f"Step2_fitting_{timestamp}.log")
logfile  = open(path_log, "w")
subprocess.run( ["conda", "run", "-n", "SwissKnife", "python", 
                 os.path.join(cfg['code_path'], 'processing_dmrs','Step2_fitting.py')] 
                + [cfg['data_path']] , env=env,stdout=logfile,stderr=logfile, check=True)
logfile.close()