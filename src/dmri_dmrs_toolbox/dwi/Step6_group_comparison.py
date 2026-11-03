"""Average subject model parameter maps already registered to atlas space."""

import fnmatch
import glob
import os
import nibabel as nib
import numpy as np
import pandas as pd
from dmri_dmrs_toolbox.misc.bids_structure import create_bids_structure
from dmri_dmrs_toolbox.misc.custom_functions_models import get_param_names_model
from dmri_dmrs_toolbox.misc.atlas_functions import prepare_atlas


def Step6_group_comparison(cfg):
    """Average subject parameter maps already registered to atlas space."""
    
    data_path = cfg['data_path']
    _, template_path = prepare_atlas(cfg['atlas'], cfg['common_folder'],'atlas', cfg)

    atlas_reference = nib.load(template_path)

    scan_list = pd.read_excel(os.path.join(data_path, cfg['scan_list_name']))
    selected = scan_list.loc[
        scan_list['study_name'].isin(cfg['subj_list'])
        & (scan_list['analyse'] == 'y'),
        ['study_name', 'group', 'sessNo'],
    ].dropna(subset=['study_name', 'group', 'sessNo']).drop_duplicates()

    configured_sessions = cfg.get('sess_list')
    if configured_sessions is None:
        sessions = sorted(selected['sessNo'].astype(int).unique())
    else:
        sessions = sorted({int(session) for session in configured_sessions})

    for sess in sessions:
        session_rows = selected[selected['sessNo'].astype(int) == sess]
        print(f'Working on session {sess}...')

        for group, group_rows in session_rows.groupby('group'):
            subjects = sorted(group_rows['study_name'].unique())
            print(f'Processing group {group}...')
            
            models_DKI=[]
            # Write list of models
            filtered_data = scan_list[
                (scan_list['phaseDir'] == 'fwd') &
                (scan_list['sessNo'] == sess) &
                (scan_list['noBval'] > 1) &
                (scan_list['acqType'] == 'PGSE') &
                (scan_list['analyse'] == 'y') 
            ]
            Delta_list = sorted(filtered_data['diffTime'].dropna().astype(int).unique())
            for Delta in Delta_list:
                models_DKI.append(f'DTI_DKI_Delta_{Delta}')
                
            for model in cfg['model_list_GM'] + cfg['model_list_WM'] +  models_DKI + ['microFA']:  
                print(f'Working with {model}...')
  
                model_clean = "DTI_DKI" if "DTI_DKI" in model else "Nexi" if "Nexi" in model else "Sandi" if "Sandi" in model else model
                patterns, _, _ = get_param_names_model(model_clean, cfg['is_alive'])
                group_path = os.path.join(
                    data_path,
                    'derivatives',
                    'group_analysis',
                    f'ses-{sess:02d}',
                    model,
                )
                os.makedirs(group_path, exist_ok=True)

                atlas_maps_by_name = {}
                for subj in subjects:
                    bids_analysis = create_bids_structure(
                        subj=subj,
                        sess=sess,
                        datatype='dwi',
                        root=data_path,
                        folderlevel='derivatives',
                        workingdir=cfg['analysis_foldername'],
                        description=model,
                    )
                    output_path = bids_analysis.get_path()

                    if model == 'microFA':
                        bids_analysis_STE = create_bids_structure(
                            subj=subj,
                            sess=sess,
                            datatype='dwi_STE',
                            root=data_path,
                            folderlevel='derivatives',
                            workingdir=cfg['analysis_foldername'],
                            description=model,
                        )        
                        output_path = bids_analysis_STE.get_path()
                        if not os.path.exists(output_path):
                            continue
                                            
                    source_path = os.path.join(output_path, 'Output_in_atlas')

                    filenames = sorted(
                        filename
                        for filename in os.listdir(source_path)
                        if any(fnmatch.fnmatch(filename, pattern) for pattern in patterns)
                    )
                    for filename in filenames:
                        print(f'For {subj}: using {filename}')
                        atlas_maps_by_name.setdefault(filename, []).append(
                            os.path.join(source_path, filename)
                        )

                for filename, atlas_maps in atlas_maps_by_name.items():
                    images = [nib.load(path) for path in atlas_maps]
                    reference_shape = atlas_reference.shape
                    reference_affine = atlas_reference.affine
                    data_arrays = []

                    for image, image_path in zip(images, atlas_maps):
                        if image.shape != reference_shape or not np.allclose(
                            image.affine, reference_affine
                        ):
                            raise ValueError(
                                f'Atlas-space map does not match the atlas template grid: {image_path}'
                            )
                        data_arrays.append(np.asarray(image.dataobj, dtype=np.float64))

                    average = np.mean(np.stack(data_arrays, axis=0), axis=0)
                    output_image = nib.Nifti1Image(
                        average.astype(np.float32),
                        reference_affine,
                        atlas_reference.header.copy(),
                    )
                    output_image.to_filename(
                        os.path.join(group_path, f'mean_{group}_{filename}')
                    )
