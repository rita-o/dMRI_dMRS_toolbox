import os
import shlex
import subprocess
import time
import numpy as np
from pathlib import Path
from joblib import Parallel, delayed
import nibabel as nib
import plotly.express as px
import plotly.graph_objects as go

from uGUIDE.config_utils import load_config_uGUIDE
from uGUIDE.estimation import estimate_microstructure


# ============================================================
# Helper functions
# ============================================================

def postprocess_NEXI_SANDIX(samples, config, adjust_parameter_limits=None):
    """
    Post-process uGUIDE samples:
    - convert U0/U1 into Di/De
    - for SANDIX, convert f/fs into neurite and soma volume fractions
    """
    prior_keys = np.array(list(config["prior"].keys()))

    idx_u0 = np.where(prior_keys == "U0")[0][0]
    idx_u1 = np.where(prior_keys == "U1")[0][0]

    u0 = samples[:, idx_u0]
    u1 = samples[:, idx_u1]

    u0 = np.clip(u0, 0, 1)
    u1 = np.clip(u1, 0, 1)

    Di_min, Di_max = config["prior_postprocessing"]["Di"]

    Di = np.sqrt((Di_max - Di_min) ** 2 * u0) + Di_min
    De = Di_min + (Di - Di_min) * u1

    out_samples = samples.copy()
    out_samples[:, idx_u0] = Di
    out_samples[:, idx_u1] = De

    if config["microstructure_model_name"].lower() == "sandix":
        idx_f = np.where(prior_keys == "f")[0][0]
        idx_fs = np.where(prior_keys == "fs")[0][0]

        f = samples[:, idx_f]
        fs = samples[:, idx_fs]

        f_neurites = f * (1 - fs)
        f_soma = f * fs

        out_samples[:, idx_f] = f_neurites
        out_samples[:, idx_fs] = f_soma

    return out_samples



def get_slice_and_masks(res, mask_map, mask_degeneracy, slice_idx, ax_plot, param_idx):
    """
    Extract a 2D slice for the chosen axis and parameter index.
    """
    if ax_plot == 0:
        res_slice = res[slice_idx, :, :, param_idx].copy()
        mask_slice = mask_map[slice_idx, :, :, param_idx]
        deg_slice = mask_degeneracy[slice_idx, :, :, param_idx].copy()
    elif ax_plot == 1:
        res_slice = res[:, slice_idx, :, param_idx].copy()
        mask_slice = mask_map[:, slice_idx, :, param_idx]
        deg_slice = mask_degeneracy[:, slice_idx, :, param_idx].copy()
    elif ax_plot == 2:
        res_slice = res[:, :, slice_idx, param_idx].copy()
        mask_slice = mask_map[:, :, slice_idx, param_idx]
        deg_slice = mask_degeneracy[:, :, slice_idx, param_idx].copy()
    else:
        raise ValueError("ax_plot must be 0, 1, or 2")

    res_slice[mask_slice == 0] = np.nan
    deg_slice[mask_slice == 0] = 0

    return res_slice, mask_slice, deg_slice


def save_map(data, dwi_img, path):
    # Publish each file after writing it so it can be opened during fitting.
    tmp_path = path.with_name(f".tmp_{path.name}")
    nib.Nifti1Image(data, affine=dwi_img.affine).to_filename(tmp_path)
    tmp_path.replace(path)
 
    
# ============================================================
# Main function
# ============================================================

def main_model_fit(main_folder, results_folder, dwi_path, mask_path,
                   model="Nexi", noise="rician", hidden_layers=[50, 30], nb_simu=700_000, nb_theta=1_000):
   
    name_extension = ""
    ax_plot = 2
    n_jobs = 7
    import torch

    if model.lower() == "nexi":
        nf_features = 14
    elif model.lower() == "sandix":
        nf_features = 22
    else:
        raise ValueError(f"Unknown model: {model}")
    
    train_size = int(nb_simu)
    
    config_file = (
        main_folder
        / "uGUIDE_config"
        / "results"
        / model
        / noise
        / f"train_size_{train_size}"
        / f"MLP_{nf_features}_features_layers_{hidden_layers[0]}_{hidden_layers[1]}"
        / f"config_{model}.pkl"
    )

    config = load_config_uGUIDE(savefile=config_file)

    # ============================================================
    # Paths
    # ============================================================
    # subject_folder = main_folder / f"{subject_id}" /  f"{session_id}" 
    # dwi_folder = subject_folder / "dwi"
    # uguide_folder = dwi_folder / f"{model}_uGUIDE"
    #results_folder = uguide_folder 
    plot_folder = results_folder / "plots"

    results_folder.mkdir(parents=True, exist_ok=True)
    plot_folder.mkdir(parents=True, exist_ok=True)

    # dwi_path = uguide_folder / "inputs" / "powderaverage_dwi.nii.gz"
    # mask_path = uguide_folder / "inputs" / "mask_dil.nii.gz"

    map_file = results_folder / f"uGUIDE_MAP_{model}{name_extension}.nii.gz"
    mask_file = results_folder / f"uGUIDE_mask_{model}{name_extension}.nii.gz"
    degeneracy_file = results_folder / f"uGUIDE_mask_degeneracy_{model}{name_extension}.nii.gz"
    uncertainty_file = results_folder / f"uGUIDE_uncertainty_{model}{name_extension}.nii.gz"
    ambiguity_file = results_folder / f"uGUIDE_ambiguity_{model}{name_extension}.nii.gz"

    print("=" * 80)
    #print(f"Subject {subject_id}, session {session_id}")
    print(f"DWI path       : {dwi_path}")
    print(f"Mask path      : {mask_path}")
    print(f"Results folder : {results_folder}")
    print("=" * 80)

    if not dwi_path.exists():
        raise FileNotFoundError(f"DWI file not found: {dwi_path}")
    if not mask_path.exists():
        raise FileNotFoundError(f"Mask file not found: {mask_path}")

    dwi_img = nib.load(dwi_path)
    dwi_data = dwi_img.get_fdata()

    mask_img = nib.load(mask_path)
    mask_data = mask_img.get_fdata()

    middle = dwi_data.shape[ax_plot] // 2

    all_map_files_exist = all(
        [
            map_file.exists(),
            mask_file.exists(),
            degeneracy_file.exists(),
            uncertainty_file.exists(),
            ambiguity_file.exists(),
        ]
    )

    # ============================================================
    # Either load existing maps or compute them
    # ============================================================
    completed_file = results_folder / f"uGUIDE_completed_{model}{name_extension}.nii.gz"
    n_params = len(config["prior_postprocessing"])
    map_shape = (*mask_data.shape, n_params)
 
    map_specs = [
        (map_file, np.float32),
        (mask_file, np.uint8),
        (degeneracy_file, np.uint8),
        (uncertainty_file, np.float32),
        (ambiguity_file, np.float32),
    ]
 
    # Select valid voxels and sort by z, then x, then y.
    i, j, k = np.nonzero(mask_data)
    valid = np.isfinite(dwi_data[i, j, k, :]).all(axis=1)
    i, j, k = i[valid], j[valid], k[valid]
    order = np.lexsort((j, i, k))
    i, j, k = i[order], j[order], k[order]
 
    # If complete map exists, get the corresponding indices and compute the indices that are missing
    if completed_file.exists():
        completed = np.array(nib.load(completed_file).dataobj, dtype=np.uint8)
        maps = [
            np.array(nib.load(path).dataobj, dtype=dtype)
            for path, dtype in map_specs
        ]
        if completed.shape != mask_data.shape or any(
            result.shape != map_shape for result in maps
        ):
            raise ValueError("Checkpoint shapes do not match the current data.")
    else:
        completed = np.zeros(mask_data.shape, dtype=np.uint8)
        maps = [np.zeros(map_shape, dtype=dtype) for _, dtype in map_specs]
    remaining = np.flatnonzero(completed[i, j, k] == 0)
    already_completed = len(i) - len(remaining)
 
    # Settings
    config["voxel_batch_size"] = 16
    config["nb_samples"] = 10_000
    n_chunks = 10
    processing_chunk_size = max(1, (len(i) + n_chunks - 1) // n_chunks)
     
    # Process chunks of data
    for start in range(0, len(remaining), processing_chunk_size):
        positions = remaining[start:start + processing_chunk_size]
        coords = (i[positions], j[positions], k[positions])
        signals = torch.as_tensor(dwi_data[coords], dtype=torch.float32)
 
        # Estimate microstructure
        estimates = estimate_microstructure(
            signals,
            config,
            postprocessing=postprocess_NEXI_SANDIX,
            plot=False,
        )
 
        for result, estimate in zip(maps, estimates):
            result[coords] = estimate.detach().cpu().numpy()
            
        # Save maps
        for result, (path, _) in zip(maps, map_specs):
            save_map(result, dwi_img, path)
 
        # Record completion only after all five maps have been saved.
        completed[coords] = 1
        save_map(completed, dwi_img, completed_file)
 
        done = already_completed + start + len(positions)
        print(f"Checkpoint saved: {done:,}/{len(i):,} voxels", flush=True)
 
    (
        param_map_nii,
        mask_map_nii,
        degeneracy_map_nii,
        uncertainty_map_nii,
        ambiguity_map_nii,
    ) = maps

    # ============================================================
    # Plot results
    # ============================================================
    nb_voxels_in_mask_map = np.count_nonzero(mask_data)

    for res, res_name in zip(
        [param_map_nii, uncertainty_map_nii, ambiguity_map_nii],
        ["map", "uncertainty", "ambiguity"],
    ):
        for p, param in enumerate(config["prior_postprocessing"].keys()):

            masked_res, _, masked_deg = get_slice_and_masks(
                res=res,
                mask_map=mask_map_nii,
                mask_degeneracy=degeneracy_map_nii,
                slice_idx=middle,
                ax_plot=ax_plot,
                param_idx=p,
            )

            valid_values = masked_res[~np.isnan(masked_res)]
            if valid_values.size == 0:
                print(f"{res_name} {param}: no valid values, skipping plot")
                continue

            print(
                f"{res_name} {param}: "
                f"min = {valid_values.min()} - max = {valid_values.max()}"
            )

            if res_name == "map":
                zmin = config["prior_postprocessing"][param][0]
                zmax = config["prior_postprocessing"][param][1]
            else:
                zmin = 0
                zmax = 50

            fig = px.imshow(
                masked_res.T,
                zmin=zmin,
                zmax=zmax,
                origin="lower",
                color_continuous_scale="viridis",
            )

            idx_deg = np.nonzero(masked_deg.T)
            fig.add_trace(
                go.Scatter(
                    x=idx_deg[1],
                    y=idx_deg[0],
                    mode="markers",
                    marker_color="red",
                    showlegend=False,
                )
            )

            fig.update_layout(
                plot_bgcolor="rgba(255, 255, 255, 0)",
                paper_bgcolor="rgba(255, 255, 255, 1)",
                font={"size": 24},
                coloraxis_colorbar_x=0.87,
            )
            fig.update_xaxes(showticklabels=False)
            fig.update_yaxes(showticklabels=False)

            plot_file = plot_folder / f"plot_{res_name}_{param}_masked_{model}_ax{ax_plot}.png"
            fig.write_image(plot_file)
            print(f"Saved plot: {plot_file}")

    # ============================================================
    # Degeneracy summary
    # ============================================================
    for p, param in enumerate(config["prior_postprocessing"].keys()):
        deg_mask_param = degeneracy_map_nii[..., p].copy()
        deg_mask_param[mask_map_nii[..., p] == 0] = 0
        nb_deg = np.count_nonzero(deg_mask_param)
        print(f"{param}: {nb_deg}/{nb_voxels_in_mask_map} degeneracies")
