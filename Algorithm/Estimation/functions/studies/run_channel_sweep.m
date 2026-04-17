function study = run_channel_sweep(cfg_base, settings)
%RUN_CHANNEL_SWEEP Compare one-satellite and all-visible multichannel modes.

counts = settings.channel_counts(:)';
n_counts = numel(counts);
n_mc = settings.n_mc;

study.channel_counts = counts;
study.metrics = cell(n_counts, n_mc);
study.mean_RMSE_h_m = nan(1, n_counts);
study.std_RMSE_h_m = nan(1, n_counts);
study.mean_NIS = nan(1, n_counts);
study.frac_converged = nan(1, n_counts);

for ci = 1:n_counts
    n_ch = counts(ci);
    fprintf('  Channels: %d (%d MC seeds)\n', n_ch, n_mc);

    rmse = nan(n_mc, 1);
    nis = nan(n_mc, 1);
    conv = nan(n_mc, 1);

    for mc = 1:n_mc
        cfg = cfg_base;
        cfg.duration_sec = settings.duration_sec;
        cfg.dt = settings.dt;
        cfg.sigma_doppler_Hz = settings.sigma_doppler_Hz;
        cfg.sigma_dopplerRate_Hzps = settings.sigma_dopplerRate_Hzps;
        cfg.seedMeas = settings.base_seed + 100*ci + mc;
        cfg.meas.useDoppler = true;
        cfg.meas.useDopplerRate = true;
        cfg.tropo.enable = false;

        if n_ch == 1
            cfg.capture.mode = 'stickyBest';
            cfg.capture.minDwell_s = 30;
            cfg.capture.switchBlackout_s = 3;
        else
            cfg.capture.mode = 'allVisible';
            cfg.maxSats = n_ch;
        end

        try
            [~, m] = run_case(cfg);
            study.metrics{ci, mc} = m;
            rmse(mc) = m.RMSE_h_m;
            nis(mc) = m.meanNIS;
            conv(mc) = m.convTime_horiz_s;
        catch ME
            warning('run_channel_sweep:runFailed', ...
                'Channels=%d seed=%d failed: %s', n_ch, mc, ME.message);
        end
    end

    study.mean_RMSE_h_m(ci) = mean(rmse, 'omitnan');
    study.std_RMSE_h_m(ci) = std(rmse, 'omitnan');
    study.mean_NIS(ci) = mean(nis, 'omitnan');
    study.frac_converged(ci) = mean(isfinite(conv) & conv < Inf, 'omitnan');
end
end
