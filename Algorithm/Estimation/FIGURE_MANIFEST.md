# Figure Manifest

This manifest links the figures used in the report PDF to the MATLAB scripts or
manual sources that produced them. It is intended for traceability; the cleaned
`main_*.m` scripts are the preferred entry points for rerunning report-facing
analyses.

## Main Report

| Report figure | Figure file | Source |
|---|---|---|
| Figure 3.1 | `Dataflow.png` | Manual work-package dataflow schematic. |
| Figure 7.1 | `Baseline_EastNorth_Errors.png` | Baseline single-satellite analysis from `legacy_reference/ANALYSIS.m` / `legacy_reference/ANALYSIS2.m`. |
| Figure 7.2 | `Baseline_ClockBias_Error.png`, `Baseline_ClockDrift_Error.png` | Baseline single-satellite clock-state analysis. |
| Figure 7.3 | `P4A_2D_Comparison.png` | Doppler-only vs Doppler+Rate comparison from Phase 4 analysis. |
| Figure 7.4 | `P4B_SmoothContourMap_SS2D_Log.png` | Noise-sensitivity sweep from `legacy_reference/ANALYSIS2.m`. |
| Figure 7.5 | `Error for differnet times PLANAR contours(4focus).png` | Time-snapshot noise-sensitivity figure from `legacy_reference/ANALYSIS2.m`. |
| Figure 8.1 | `P8 Skyplot.png`, `P8 2D Horiz Error +-1sigma.png` | Sequential single-channel study from `legacy_reference/PHASE8_FairComparison.m`. |
| Figure 8.2 | `P6_2A_MultiChannel_Overlay.png`, `P6_2B_SSError_vs_Channels_Clean.png` | Simultaneous multi-channel study from `legacy_reference/ANALYSIS2.m`; cleaned rerun entry point: `main_multisat_study.m`. |
| Figure 8.3 | `P6_2E.png`, `Filled ellipses.png` | Spatial/mapping visualisation from `legacy_reference/ANALYSIS2.m`. |

## Appendices

| Report figure | Figure file | Source |
|---|---|---|
| Figure D.1 | `Progressive Complexity.png` | Phase 0 validation workflow, `main_phase0_validation.m` / `legacy_reference/MAIN_phase0_validation.m`. |
| Figure D.2 | `Measurements_rerpresentation.png` | Phase 0 measurement-stream validation. |
| Figure D.3 | `NIS Diagnostics.png` | NIS diagnostics from the validation workflow. |
| Figure D.4 | `Doppler innovation.png` | Doppler innovation diagnostics. |
| Figure D.5 | `Doppler-rate innovation.png` | Doppler-rate innovation diagnostics. |
| Figure D.6 | `B3_Rinflate_sweep.png` | Measurement-inflation sweep from Phase 0 validation. |
| Figure D.7 | `B4_P0_sweep.png` | Initial-covariance sweep from Phase 0 validation. |
| Figure D.8 | `B5_process_noise_sweep.png` | Process-noise sweep from Phase 0 validation. |
| Figure E.1 | `Baseline_2D_Horizontal.png` | Doppler-only 2D baseline result. |
| Figure E.2 | `Convergence Topographic (converging time).png` | Focus-region convergence-time map from the noise sweep. |
| Figure E.3 | `BestCaseConvergenceSingleSatellie.png` | Best-case single-satellite convergence run. |
| Figure E.4 | `Bestcaseclock.png` | Best-case clock-state result. |
| Figure F.1 | `P2A_2D_Horizontal_Tropo.png` | Tropospheric sensitivity study. |
| Figure F.2 | `P2D_ClkBias_Tropo.png`, `P2E_ClkDrift_Tropo.png` | Tropospheric clock-state sensitivity. |
| Figure G.1 | `APP_OW1_skyplot.png` | OneWeb geometry analysis from `legacy_reference/ANALYSIS_GEOM.m`. |
| Figure G.2 | `APP_OW4_visibility.png` | Elevation-only visibility analysis from OneWeb geometry scripts. |
| Figure G.3 | `APP_OW5_azimuth_diversity.png` | Azimuth-diversity analysis from OneWeb geometry scripts. |
| Figure G.4 | `APP_OW6_doppler_range.png` | Doppler-excursion analysis using OneWeb geometry outputs. |
| Figure G.5 | `Screenshot 2026-03-23 233204.png` | External ECC Report 271 beam-pattern screenshot used as a cited source figure. |
| Figure G.6 | `Screenshot 2026-03-23 233306.png` | External ECC Report 271 3 dB contour screenshot used as a cited source figure. |

## Notes

- Some report figures were generated during exploratory analysis and then
  copied into the LaTeX `Figures/` folder. The corresponding exploratory
  scripts are retained in `legacy_reference` for provenance.
- The cleaned `functions/` tree is toolbox-light: it avoids Statistics Toolbox
  and Signal Processing Toolbox dependencies. Parallel loops can still be
  enabled where the original scripts use them.
