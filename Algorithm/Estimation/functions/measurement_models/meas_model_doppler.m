function [zpred, H] = meas_model_doppler(x, rNom, R_n2e, ps, vs, k, s, cfg, meas, iCB, nx)

if cfg.meas.enable_lighttime
    zgeom   = compute_rhodot_lighttime(rNom, cfg.rx.vel_ecef(:), ps, vs, cfg.c);
    H_ecef6 = jac_rhodot_lighttime(rNom, cfg.rx.vel_ecef(:), ps, vs, cfg.c);
else
    zgeom   = compute_rhodot(rNom, cfg.rx.vel_ecef(:), ps, vs);
    H_ecef6 = jac_rhodot(rNom, cfg.rx.vel_ecef(:), ps, vs);
end

H_pos_ecef = H_ecef6(1,1:3);
H_pos_enu  = H_pos_ecef * R_n2e;

H      = zeros(1,nx);
H(1:3) = H_pos_enu;
H(iCB) = 1;

zpred = zgeom + meas.tropoRateModel(k,s) + x(iCB);
end
