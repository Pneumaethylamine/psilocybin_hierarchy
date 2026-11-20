# Transient compression of brain hierarchy during the psilocybin-induced psychedelic state in humans

Feature-extraction and statistics code for:

> McCulloch DE-W, Olsen AS, Larsen K, Madsen MK, Sanz Perl Y, Ozenne B, Knudsen GM, Deco G, Fisher PM, *et al.*
> Transient compression of brain hierarchy during the psilocybin-induced psychedelic state in humans. (submitted)

The input is parcellated resting-state fMRI time series: 200 Schaefer cortical regions plus 32 Tian S3
subcortical regions, already preprocessed in SPM12 and denoised in CONN. Preprocessing code is not included,
and participant-level data are not included. Data are available on request through the Cimbi database.

Run the scripts in step order.

## Feature extraction (MATLAB)

| Script | What it does |
|---|---|
| `matlab/step1_extract_insideout.m` | Chooses the lag τ for each scanner from the autocorrelation decay. Computes the time-shifted correlation and mutual-information matrices, global and regional INSIDEOUT, and regional IN/OUT for every scan. |
| `matlab/compute_insideout.m` | Core function: time-shifted MI matrix and INSIDEOUT. |
| `matlab/calculate_autocorrelations.m`, `matlab/plot_autocorrelations.m` | Regional autocorrelation curves used to choose τ. |
| `matlab/step2_modularity.m` | Modularity of the static FC and time-shifted MI matrices for each scan, with null-model distributions. |

Requires the [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) (`community_louvain`, `null_model_und_sign`, `null_model_dir_sign`).

## Statistics (R)

| Script | What it does |
|---|---|
| `R/step3_psilocin_mmrm_models.R` | MMRM models of each metric against plasma psilocin level (global and the 232 regions), hierarchy-position analyses, transmodal vs unimodal permutation test, time-window comparisons and figures. |
| `R/step4_hierarchy_bootstrap.R` | Participant-level bootstrap of the correlation between baseline hierarchy position and regional psilocin effects. |
| `R/step5_motion_scanner_timepoint.R` | Head-motion (FD) analyses, PPL × scanner interaction models, and paired t-tests of each time window against pre-drug. |

Model: `mmrm(metric ~ PPL + FD + Age + Sex + Scanner + us(scan | participant))`. Plasma psilocin is in mg/L.
R packages: `mmrm`, `dplyr`, `tidyr`, `ggplot2`, `readr`, `readxl`, `writexl`, `scales`, `LMMstar`, `MuMIn`.

The scripts also include analyses of global functional connectivity (GCOR).
File locations are read from environment variables at the top of each script (`UPF_PROJECT_DIR`, `NP2_DATA_DIR`,
`BCT_DIR`, `ROI_DIR`).

## Licence

See `LICENSE`.
