function H = jac_rhodot(r_r, ~, r_s, v_s)
dr  = r_s - r_r;
dv  = v_s;
rho = norm(dr);
drhod_dp = (-dv/rho + (dv.'*dr)*dr/(rho^3)).';
drhod_dv = (-dr/rho).';
H = [drhod_dp, drhod_dv];
end
