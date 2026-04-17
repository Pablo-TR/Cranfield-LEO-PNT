# LEO Doppler EKF Report Code

## Folder Structure

```text
report code/
    InputData/
        oneweb.tle
    main_baseline.m
    main_phase0_validation.m
    main_noise_sensitivity.m
    main_multisat_study.m
    functions/
        config/
            cfg_report_baseline.m
        metrics/
            acf_manual.m
            ekf_extract_metrics.m
            extract_metrics.m
            find_convergence_time.m
            print_metric_summary.m
        propagation/
            read_tle_file_fast.m
            propagate_tle_ecef_fast.m
            propagate_j2.m
            eci2ecef.m
        constellation/
            build_constellation_fast.m
        geometry/
            compute_*.m
            jac_*.m
            geodetic2ecef_local.m
            wrapTo360.m
        measurements/
            generate_measurements_fastrealistic.m
        measurement_models/
            meas_model_*.m
        filter/
            scalar_iekf_update.m
            vector_iekf_update.m
        atmosphere/
            tropo_delay_simple.m
            ekf_tropo_ztd_from_met.m
        plotting/
            ekf_plot_*.m
            plot_baseline_report.m
            plot_noise_sensitivity.m
            plot_multisat_summary.m
        solver/
            real_strtg_sequential_EKF_fastrealistic.m
            real_strtg_sequential_EKF_geom.m
            ekf_tropo_dispatch.m
            leo_iekf.m
            leo_geom_iekf.m
        studies/
            run_case.m
            run_phase0_validation.m
            define_noise_sweeps.m
            run_sweep1d.m
            run_sweep2d.m
            run_channel_sweep.m
        utils/
            setup_report_paths.m
            report_root.m
            ensure_dir.m
            set_field.m
            safe_mean.m
            safe_std.m
        legacy_reference/
        ANALYSIS.m
        ANALYSIS2.m
        ANALYSIS_GEOM.m
        GEN_DOPPLER_CURVES.m
        MAIN_phase0_validation.m
        MAIN_sensitivity.m
        PHASE8_FairComparison.m
    FIGURE_MANIFEST.md
```

## Main Scripts

- `main_baseline.m` runs one adopted baseline case and saves a compact
  diagnostic figure and `.mat` file.
- `main_phase0_validation.m` runs the full B1-B6 validation workflow preserved
  from `MAIN_phase0_validation.m`.
- `main_noise_sensitivity.m` runs the cleaned 1D and 2D parameter sensitivity
  sweeps from the original sensitivity scripts.
- `main_multisat_study.m` runs a cleaned multichannel study based on the
  `ANALYSIS2` multichannel logic.

Run scripts from inside `report code` or open them directly in MATLAB. Each
script adds the local `functions/` tree to the path and uses the copied
`InputData/oneweb.tle`.

## Execution Order

1. `main_baseline.m`
2. `main_phase0_validation.m`
3. `main_noise_sensitivity.m`
4. `main_multisat_study.m`

The sensitivity and multichannel scripts run many solver calls. Reduce
`settings.n_repeat`, `settings.n_mc`, or sweep arrays in the top-level scripts
for quick smoke tests.

## Refactoring Summary

- Centralised repeated configuration into `cfg_report_baseline.m` with profiles
  for `phase0`, `sensitivity`, `baseline300`, and `multisat`.
- Added short-name solver wrappers `leo_iekf.m` and `leo_geom_iekf.m`.
- Separated orbit/TLE propagation from the solvers into
  `functions/propagation`.
- Moved shared constellation, measurement, measurement-model, IEKF-update,
  geometry, atmosphere, plotting, and metric helpers out of the solver files.
- Added `run_case.m`, `run_sweep1d.m`, `run_sweep2d.m`, and
  `run_channel_sweep.m` so Monte Carlo and sweep logic is no longer embedded
  directly in every top-level script.
- Added a defensive `extract_metrics.m` wrapper for cases where Doppler-rate
  innovations are intentionally absent.
- Replaced Signal Processing Toolbox `xcorr` use with the toolbox-free
  `acf_manual.m`. No Statistics Toolbox functions are required by the cleaned
  report-facing code.
- Added `FIGURE_MANIFEST.md` to trace report figures to their generating
  scripts or manual sources.

## Intentionally Unchanged

- The core IEKF equations and numerical settings.
- The geometry/handover policy logic.
- Existing Phase 0 plotting and console text inside the preserved validation
  script.

## Known Issues and Engineering Notes

- The original workspace contains several baseline configurations with
  different durations, Doppler-rate noise values, troposphere settings, and
  clock process-noise values. These are now explicit profiles in
  `cfg_report_baseline.m` rather than hidden local functions.
- `ANALYSIS2.m` contains multiple later-phase experiments and patch notes in a
  single long script. It is preserved in `legacy_reference`; only the
  multichannel workflow needed for the report-facing script was extracted.
- Legacy scripts may contain old comments and exploratory branches. They are
  retained only for provenance; use the `main_*.m` scripts for report-facing
  execution.
- The full Phase 0 workflow is computationally heavy because B2-B6 use Monte
  Carlo runs. It is not run automatically by the other scripts.

