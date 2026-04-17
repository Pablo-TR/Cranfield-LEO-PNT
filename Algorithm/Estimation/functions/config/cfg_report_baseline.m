function cfg = cfg_report_baseline(profile)
%CFG_REPORT_BASELINE Shared configuration builder for report scripts.
%
% Profiles:
%   phase0      adopted validation baseline used by MAIN_phase0_validation
%   sensitivity compact sensitivity profile from MAIN_sensitivity
%   baseline300 single-run demonstration values from MAIN_baseline
%   multisat    multichannel values from ANALYSIS2 / Phase 8

if nargin < 1 || isempty(profile)
    profile = 'phase0';
end

root_dir = report_root();
cfg.tleFile = fullfile(root_dir, 'InputData', 'oneweb.tle');
cfg.duration_sec = 100;
cfg.dt = 0.01;
cfg.startUTC = [];
cfg.stopUTC = [];
cfg.startFromTLE = 'first';
cfg.maxTLEAgeDays = 7;

cfg.rx.lat_deg = 52.069;
cfg.rx.lon_deg = -0.609;
cfg.rx.alt_m = 100;
cfg.rx.vel_ecef = [0;0;0];
cfg.rx.acc_ecef = [0;0;0];

cfg.fc_Hz = 11.325e9;
cfg.c = 299792458;
cfg.lambda = cfg.c / cfg.fc_Hz;

cfg.minElevDeg = 25;
cfg.coarseStep = 10;
cfg.minMaxElev = 30;
cfg.maxSats = 12;
cfg.constellationSelectionWindow_s = 60;

cfg.useSGP4 = false;
cfg.allowJ2Fallback = true;

cfg.init.errENU_m = [400; 400; 0];
cfg.init.sigENU_m = [500; 500; 10];
cfg.P0_clkBias = (0.2)^2;
cfg.P0_clkDrift = (5e-7)^2;
cfg.sigma_pos_rw = 1e-4;

cfg.clk.truth.initBias_mps = 0.10;
cfg.clk.truth.initDrift_mps2 = 5e-5;
cfg.clk.truth.sigmaBias_rw_mps_sqrt = 1e-2;
cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 3e-4;
cfg.clk.filt.sigmaBias_rw_mps_sqrt = 1e-2;
cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = 5e-4;
cfg.clk.filter_init_bias_mps = 0.0;
cfg.clk.filter_init_drift_mps2 = 0.0;

cfg.tropo.enable = true;
cfg.tropo.ztd_truth_m = 2.5;
cfg.tropo.ztd_model_m = 2.5;
cfg.tropo.el_floor_deg = 5;
cfg.tropo.map_epsilon_deg = 0.5;
cfg.tropo.model_truth = 'simple';
cfg.tropo.model_filter = 'simple';
cfg.tropo.lat_deg = cfg.rx.lat_deg;
cfg.tropo.alt_m = cfg.rx.alt_m;
cfg.tropo.doy = 180;
cfg.tropo.P0_hPa = 1013.25;
cfg.tropo.T0_K = 288.15;
cfg.tropo.e0_hPa = 11.0;

cfg.sigma_doppler_Hz = 50.0;
cfg.sigma_dopplerRate_Hzps = 1.0;
cfg.meas.enable_lighttime = false;
cfg.meas.colored_tau_s = 0.5;
cfg.meas.sigma_col0_mps = 0.03;
cfg.meas.Rinflate = 2.0;
cfg.meas.useDoppler = true;
cfg.meas.useDopplerRate = true;
cfg.meas.fdJacStep_m = 0.10;

cfg.filter.useNISGate = true;
cfg.filter.nisGate = 16.0;
cfg.filter.maxIEKFIter = 10;
cfg.filter.iekfTol = 1e-5;
cfg.filter.maxStep_m = 25;
cfg.filter.startupInflateEpochs = 50;
cfg.filter.startupRmult = 25;

cfg.height.enable = true;
cfg.height.sigma_m = 10.0;
cfg.height.z0_m = 0.0;

cfg.capture.mode = 'stickyBest';
cfg.capture.minDwell_s = 30;
cfg.capture.switchBlackout_s = 3;
cfg.capture.elev_margin_deg = 20;
cfg.capture.useServiceAngleFilter = false;
cfg.capture.maxServiceOffNadirDeg = 45;
cfg.capture.switchThreshold_lowElev = 10;
cfg.capture.switchThreshold_highElev = 15;
cfg.capture.highElevThreshold = 35;

cfg.analysis.conv.refWindow_s = 15;
cfg.analysis.conv.settleWindow_s = 10;
cfg.analysis.conv.relTol = 0.15;
cfg.analysis.conv.absTol_m = 5.0;
cfg.analysis.conv.maxCoeffVar = 0.25;
cfg.analysis.conv.maxSlope_mps = 1.0;
cfg.analysis.conv.errFloor_m = 5.0;

cfg.seedMeas = 42;
cfg.plot.enable = false;
cfg.verbose = false;

switch lower(profile)
    case 'phase0'
        % Adopted validation baseline above.
    case 'baseline300'
        cfg.duration_sec = 300;
        cfg.rx.lat_deg = 52.06869;
        cfg.rx.lon_deg = -0.60884;
        cfg.sigma_dopplerRate_Hzps = 3.0;
        cfg.clk.truth.sigmaBias_rw_mps_sqrt = 5e-3 * 2;
        cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4 * 2;
    case 'sensitivity'
        cfg.duration_sec = 5;
        cfg.rx.lat_deg = 52.06869;
        cfg.rx.lon_deg = -0.60884;
        cfg.init.errENU_m = [500;500;0];
        cfg.init.sigENU_m = [550;550;10];
        cfg.P0_clkBias = (0.1)^2;
        cfg.P0_clkDrift = (5e-5)^2;
        cfg.clk.truth.sigmaBias_rw_mps_sqrt = 5e-3;
        cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
        cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = 1e-4;
        cfg.tropo.ztd_model_m = 2.3;
        cfg.sigma_doppler_Hz = 5.0;
        cfg.sigma_dopplerRate_Hzps = 0.05;
        cfg.meas.enable_lighttime = true;
        cfg.meas.colored_tau_s = 30;
        cfg.meas.Rinflate = 1.0;
        cfg.filter.startupInflateEpochs = 20;
    case 'multisat'
        cfg.duration_sec = 100;
        cfg.dt = 0.01;
        cfg.rx.lat_deg = 52.06869;
        cfg.rx.lon_deg = -0.60884;
        cfg.tropo.enable = false;
        cfg.clk.truth.sigmaBias_rw_mps_sqrt = 5e-3;
        cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
        cfg.P0_clkBias = (0.15)^2;
        cfg.P0_clkDrift = (2.5e-7 * 2)^2;
    otherwise
        error('Unknown report config profile: %s', profile);
end
end
