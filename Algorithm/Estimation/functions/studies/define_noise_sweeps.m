function sweeps = define_noise_sweeps()
%DEFINE_NOISE_SWEEPS Parameter sweeps used by main_noise_sensitivity.

i = 0;

i = i + 1;
sweeps(i).param_path = 'sigma_pos_rw';
sweeps(i).values = logspace(-5, -2, 12);
sweeps(i).label = 'sigma_pos_rw';
sweeps(i).units = 'm/sqrt(s)';
sweeps(i).tag = 'sigma_pos_rw';

i = i + 1;
sweeps(i).param_path = 'sigma_doppler_Hz';
sweeps(i).values = logspace(0, 2, 12);
sweeps(i).label = 'sigma_doppler_Hz';
sweeps(i).units = 'Hz';
sweeps(i).tag = 'sigma_doppler_Hz';

i = i + 1;
sweeps(i).param_path = 'sigma_dopplerRate_Hzps';
sweeps(i).values = logspace(-2, 1, 12);
sweeps(i).label = 'sigma_dopplerRate_Hzps';
sweeps(i).units = 'Hz/s';
sweeps(i).tag = 'sigma_dopplerRate_Hzps';

i = i + 1;
sweeps(i).param_path = 'meas.Rinflate';
sweeps(i).values = [0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12];
sweeps(i).label = 'Rinflate';
sweeps(i).units = '-';
sweeps(i).tag = 'Rinflate';

i = i + 1;
sweeps(i).param_path = 'height.sigma_m';
sweeps(i).values = [1, 2, 5, 10, 20, 50, 100, 500, 1000];
sweeps(i).label = 'height_sigma_m';
sweeps(i).units = 'm';
sweeps(i).tag = 'height_sigma';

i = i + 1;
sweeps(i).param_path = 'P0_clkBias';
sweeps(i).values = logspace(-4, 0, 10).^2;
sweeps(i).label = 'sqrt_P0_clkBias';
sweeps(i).units = 'm/s';
sweeps(i).tag = 'P0_clkBias';

i = i + 1;
sweeps(i).param_path = 'P0_clkDrift';
sweeps(i).values = logspace(-6, -2, 10).^2;
sweeps(i).label = 'sqrt_P0_clkDrift';
sweeps(i).units = 'm/s^2';
sweeps(i).tag = 'P0_clkDrift';

i = i + 1;
sweeps(i).param_path = 'clk.filt.sigmaBias_rw_mps_sqrt';
sweeps(i).values = logspace(-4, -1, 10);
sweeps(i).label = 'clock_bias_rw';
sweeps(i).units = 'm/s/sqrt(s)';
sweeps(i).tag = 'clk_bias_rw';

i = i + 1;
sweeps(i).param_path = 'clk.filt.sigmaDrift_rw_mps2_sqrt';
sweeps(i).values = logspace(-6, -2, 10);
sweeps(i).label = 'clock_drift_rw';
sweeps(i).units = 'm/s^2/sqrt(s)';
sweeps(i).tag = 'clk_drift_rw';
end
