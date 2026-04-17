function [meas, aux] = appendix_measurement_generator( ...
    sat, vis, recTrue, recVelTrue, cfg)
%APPENDIX_MEASUREMENT_GENERATOR Synthetic Doppler measurement generation.
%
% This report listing summarises the measurement generator used for the
% Monte Carlo studies. It produces Doppler and Doppler-rate observations,
% including truth clock states, tropospheric model terms, elevation-weighted
% noise, and a coloured Doppler residual.

[K,Nsat] = size(vis);

meas.dopp.true  = nan(K,Nsat);
meas.dopp.obs   = nan(K,Nsat);
meas.dopp.sigma = nan(K,Nsat);

meas.rate.true  = nan(K,Nsat);
meas.rate.obs   = nan(K,Nsat);
meas.rate.sigma = nan(K,Nsat);

if cfg.tropo.enable
    tropoDelayTruth = ekf_tropo_dispatch( ...
        cfg.tropo.model_truth, sat.elev, cfg.tropo.ztd_truth_m, cfg.tropo);
    tropoDelayModel = ekf_tropo_dispatch( ...
        cfg.tropo.model_filter, sat.elev, cfg.tropo.ztd_model_m, cfg.tropo);

    tropoRateTruth = deriv1_cols(tropoDelayTruth.', cfg.dt).';
    tropoRateModel = deriv1_cols(tropoDelayModel.', cfg.dt).';
    tropoAccTruth  = deriv1_cols(tropoRateTruth.', cfg.dt).';
    tropoAccModel  = deriv1_cols(tropoRateModel.', cfg.dt).';
else
    tropoRateTruth = zeros(K,Nsat);
    tropoRateModel = zeros(K,Nsat);
    tropoAccTruth  = zeros(K,Nsat);
    tropoAccModel  = zeros(K,Nsat);
end

clkBiasTrue  = zeros(K,1);
clkDriftTrue = zeros(K,1);

clkBiasTrue(1)  = cfg.clk.truth.initBias_mps;
clkDriftTrue(1) = cfg.clk.truth.initDrift_mps2;

for k = 2:K
    clkDriftTrue(k) = clkDriftTrue(k-1) + ...
        cfg.clk.truth.sigmaDrift_rw_mps2_sqrt * sqrt(cfg.dt) * randn();

    clkBiasTrue(k) = clkBiasTrue(k-1) + ...
        clkDriftTrue(k-1) * cfg.dt + ...
        cfg.clk.truth.sigmaBias_rw_mps_sqrt * sqrt(cfg.dt) * randn();
end

phi = exp(-cfg.dt / cfg.meas.colored_tau_s);
col = zeros(K,Nsat);
for s = 1:Nsat
    rng(cfg.seedMeas + s * 100);
    for k = 2:K
        sig_col = cfg.meas.sigma_col0_mps / sind(max(sat.elev(k,s),10));
        col(k,s) = phi*col(k-1,s) + sqrt(1-phi^2)*sig_col*randn();
    end
end

for s = 1:Nsat
    rng(cfg.seedMeas + s * 200);
    for k = find(vis(:,s)).'
        ps = sat.pos(:,1,s,k);
        vs = sat.vel(:,1,s,k);
        as = sat.acc(:,1,s,k);

        if cfg.meas.enable_lighttime
            zgeomD = compute_rhodot_lighttime(recTrue, recVelTrue, ps, vs, cfg.c);
            zgeomR = compute_rhoddot_lighttime( ...
                recTrue, recVelTrue, cfg.rx.acc_ecef(:), ps, vs, as, cfg.c);
        else
            zgeomD = compute_rhodot(recTrue, recVelTrue, ps, vs);
            zgeomR = compute_rhoddot( ...
                recTrue, recVelTrue, cfg.rx.acc_ecef(:), ps, vs, as);
        end

        el_deg = max(sat.elev(k,s), cfg.tropo.el_floor_deg);
        sigmaD = cfg.meas.sigma_rhodot0_mps   / sind(el_deg);
        sigmaR = cfg.meas.sigma_rhoddot0_mps2 / sind(el_deg);

        meas.dopp.true(k,s) = zgeomD + tropoRateTruth(k,s) + clkBiasTrue(k);
        meas.rate.true(k,s) = zgeomR + tropoAccTruth(k,s)  + clkDriftTrue(k);

        meas.dopp.obs(k,s) = meas.dopp.true(k,s) + sigmaD*randn() + col(k,s);
        meas.rate.obs(k,s) = meas.rate.true(k,s) + sigmaR*randn();

        meas.dopp.sigma(k,s) = cfg.meas.Rinflate * ...
            sqrt(sigmaD^2 + (0.5*cfg.meas.sigma_col0_mps)^2);
        meas.rate.sigma(k,s) = cfg.meas.Rinflate * sigmaR;
    end
end

meas.tropoRateModel = tropoRateModel;
meas.tropoAccModel  = tropoAccModel;

aux.clkBiasTrue    = clkBiasTrue;
aux.clkDriftTrue   = clkDriftTrue;
aux.tropoRateModel = tropoRateModel;
aux.tropoAccModel  = tropoAccModel;
end
