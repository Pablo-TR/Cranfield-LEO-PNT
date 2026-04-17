function [zpred, H] = appendix_doppler_rate_measurement_model( ...
    x, rNom, R_n2e, ps, vs, as, k, s, cfg, meas, iCB, iCD, nx)
%APPENDIX_DOPPLER_RATE_MEASUREMENT_MODEL Joint Doppler and rate model.
%
% The first row is the Doppler/range-rate model and is directly sensitive
% to receiver clock bias. The second row is the Doppler-rate/range-
% acceleration model and is directly sensitive to receiver clock drift.

if cfg.meas.enable_lighttime
    zgeom_D   = compute_rhodot_lighttime(rNom, cfg.rx.vel_ecef(:), ps, vs, cfg.c);
    H_ecef6_D = jac_rhodot_lighttime(rNom, cfg.rx.vel_ecef(:), ps, vs, cfg.c);
else
    zgeom_D   = compute_rhodot(rNom, cfg.rx.vel_ecef(:), ps, vs);
    H_ecef6_D = jac_rhodot(rNom, cfg.rx.vel_ecef(:), ps, vs);
end

H_D      = zeros(1,nx);
H_D(1:3) = H_ecef6_D(1,1:3) * R_n2e;
H_D(iCB) = 1;

zpred_D = zgeom_D + meas.tropoRateModel(k,s) + x(iCB);

if cfg.meas.enable_lighttime
    zgeom_R = compute_rhoddot_lighttime( ...
        rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ps, vs, as, cfg.c);
else
    zgeom_R = compute_rhoddot( ...
        rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ps, vs, as);
end

H_pos_ecef_R = jac_rhoddot_pos_num( ...
    rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ps, vs, as, cfg);

H_R      = zeros(1,nx);
H_R(1:3) = H_pos_ecef_R * R_n2e;
H_R(iCD) = 1;

zpred_R = zgeom_R + meas.tropoAccModel(k,s) + x(iCD);

zpred = [zpred_D; zpred_R];
H     = [H_D; H_R];
end
