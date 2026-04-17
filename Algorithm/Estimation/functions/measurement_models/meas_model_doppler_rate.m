function [zpred, H] = meas_model_doppler_rate(x, rNom, R_n2e, ps, vs, as, k, s, cfg, meas, iCD, nx)

if cfg.meas.enable_lighttime
    zgeom = compute_rhoddot_lighttime(rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ...
        ps, vs, as, cfg.c);
else
    zgeom = compute_rhoddot(rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ps, vs, as);
end

H_pos_ecef = jac_rhoddot_pos_num(rNom, cfg.rx.vel_ecef(:), cfg.rx.acc_ecef(:), ps, vs, as, cfg);
H_pos_enu  = H_pos_ecef * R_n2e;

H      = zeros(1,nx);
H(1:3) = H_pos_enu;
H(iCD) = 1;

zpred = zgeom + meas.tropoAccModel(k,s) + x(iCD);
end
