"""Find candidate organoid tissue in a 3D b0 NIfTI.

Usage:
    python segment_organoids_b0.py b0.nii.gz organoids_mask.nii.gz

Dependencies: numpy, scipy, nibabel, matplotlib.
Input: a 3D b0 with zero-valued voxels outside the flask, or provide
--flask-mask with a binary flask-interior mask on the same voxel grid.
Output: a binary NIfTI on the input grid and a PNG overlay.
Review the overlay before using the mask for quantitative analysis.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy import ndimage as ndi
import pandas as pd
import os
import math 
from dmri_dmrs_toolbox.misc.bids_structure import create_bids_structure
from dmri_dmrs_toolbox.misc.custom_functions import (create_directory, copy_file)
                                                     
def Step_segment_organoids_b0(cfg):
    """Segment the selected subject/session pairs from the preprocessing scan list.

    Optional cfg keys:
        organoid_polarity: "dark" (default), "bright", or "both".
        organoid_polarity_by_scan: {(subject, session_number): "bright", ...}.
        organoid_contrast: fractional threshold, default 0.15.
        redo_organoid_masks: set True to replace existing masks, default False.

    Returns lists of written, existing, missing and failed outputs.
    """
    data_path   = cfg['data_path']     
    scan_list   = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))
    topupon     = cfg['do_topup']
    anat_format = cfg['anat_format']
    rim_mm = 0.6
    background_mm = 0.50
    smooth_mm = 0.10
    threshold = 0.10
    min_volume_mm3 = 0.03
    max_volume_mm3 = np.inf  
    void_fraction = 0.1  
    void_radius_mm = 0.15
    
    def masked_smooth(data, mask, sigma):
        weighted = ndi.gaussian_filter(np.where(mask, data, 0.0), sigma)
        weights = ndi.gaussian_filter(mask.astype(float), sigma)
        return np.divide(weighted, weights, out=np.zeros_like(weighted), where=weights > 1e-8)
    
    ######## SUBJECT-WISE OPERATIONS ########
    for subj in cfg['subj_list']:
    
        print('Preprocessing ' + subj + '...')
    
        # Extract data for subject
        subj_data      = scan_list[(scan_list['study_name'] == subj)].reset_index(drop=True)
        subj_data      = subj_data[subj_data['analyse'] == 'y'].reset_index(drop=True)

        # List of available acquisition sessions
        all_sessions = sorted(
            int(x) for x in subj_data["sessNo"].unique() if not math.isnan(x)
        )
        
        # Use all sessions unless specific ones are requested
        if cfg["sess_list"] is None:
            sess_list = all_sessions
        else:
            sess_list = [s for s in cfg["sess_list"] if s in all_sessions]
    
        ######## SESSION-WISE OPERATIONS ########
        for sess in sess_list:
            
            subj_data_sess = subj_data[subj_data["sessNo"] == sess]
            
            print('Working on session ' + str(sess) + '...')
            
            session_id = f"ses-{sess:02d}"
            
            output_folder = os.path.join(data_path, "derivatives", "manual_masks", f"{subj}", f"{session_id}")
            create_directory(output_folder)
            output_path = os.path.join(output_folder,f"{subj}_{session_id}_mask_organoids_b0.nii.gz")

            bids_strc = create_bids_structure(
                subj=subj, sess=sess, datatype="dwi", root=str(data_path),
                folderlevel="derivatives", workingdir=cfg["prep_foldername"],
            )
            bids_strc.set_param(description="allDelta-allb")
            b0 = bids_strc.get_path("b0_dn_gc_ec_avg_bc_brain.nii.gz")


            img = nib.load(b0)
            data = img.get_fdata()
            if data.ndim != 3 or not np.isfinite(data).all():
                print("Supply a finite 3D b0 volume.")
            spatial_unit = img.header.get_xyzt_units()[0]
            if spatial_unit not in ("mm", "unknown"):
                print("The script expects voxel dimensions in millimetres.")
            spacing = np.asarray(img.header.get_zooms()[:3])
           

            flask = np.zeros(data.shape, dtype=bool)
            for z in range(data.shape[2]):
                filled = ndi.binary_fill_holes(data[:, :, z] != 0)
                components, count = ndi.label(filled)
                if count:
                    sizes = np.bincount(components.ravel())
                    flask[:, :, z] = components == (1 + np.argmax(sizes[1:]))
            if not np.any(flask & (data > 0)):
                print("The flask mask contains no positive b0 signal.")
            
            interior = np.zeros_like(flask)
            for z in range(data.shape[2]):
                padded = np.pad(flask[:, :, z], 1)
                distance = ndi.distance_transform_edt(padded, sampling=spacing[:2])[1:-1, 1:-1]
                interior[:, :, z] = flask[:, :, z] & (distance > rim_mm)
            if not interior.any():
                print("No voxels remain after excluding the rim.")
                
            n_top = 5
            occupied = np.flatnonzero(flask.any(axis=(0, 1)))
     
            if n_top:
                interior[:, :, occupied[-n_top:]] = False
         
            baseline = masked_smooth(data, flask, background_mm / spacing)
            signal = masked_smooth(data, flask, smooth_mm / spacing)
            floor = 0.1 * np.median(data[flask & (data > 0)])
            contrast_map = (signal - baseline) / np.maximum(baseline, floor)
            
            labels = np.zeros(data.shape, dtype=np.uint8)
            min_voxels = max(1, int(np.ceil(min_volume_mm3 / np.prod(spacing))))
            polarities = {"dark": ((1, -1),), "bright": ((2, 1),),
                          "both": ((1, -1), (2, 1))}
            reference_signal = np.median(data[flask & (data > 0)])
            void_core = flask & (data < void_fraction * reference_signal)
            
            void_neighbourhood = void_core.copy()
            if void_core.any():
                distance_to_void = ndi.distance_transform_edt(
                    ~void_core, sampling=spacing
                )
                void_neighbourhood = distance_to_void <= void_radius_mm
            
            void_neighbourhood = void_core.copy()
            if void_core.any():
                distance_to_void = ndi.distance_transform_edt(
                    ~void_core, sampling=spacing
                )
                void_neighbourhood = distance_to_void <= 0.30  # mmr


            for label, sign in polarities["dark"]:
                candidates = (
                        interior
                        & (sign * contrast_map >= threshold)
                        & ~void_neighbourhood
                    )
                components, _ = ndi.label(candidates)
                sizes = np.bincount(components.ravel())
                component_volumes_mm3 = sizes * np.prod(spacing)
                keep = (
                    (sizes >= min_voxels)
                    & (component_volumes_mm3 <= max_volume_mm3)
                )
                keep[0] = False
                labels[keep[components]] = label
            

            header = img.header.copy()
            header.set_data_dtype(np.uint8)
            result = nib.Nifti1Image((labels > 0).astype(np.uint8), img.affine, header)
            result.set_qform(img.get_qform(), int(img.header["qform_code"]))
            result.set_sform(img.get_sform(), int(img.header["sform_code"]))
            
            nib.save(result, output_path)
            filename = os.path.basename(b0)
            copy_file([b0],  [os.path.join(output_folder,filename)])
            print(f"Mask saved in: {output_path}")
                           
            



