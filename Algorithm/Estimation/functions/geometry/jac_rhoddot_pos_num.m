function Hpos = jac_rhoddot_pos_num(r_r, v_r, a_r, r_s, v_s, a_s, cfg)
d = cfg.meas.fdJacStep_m;
Hpos = zeros(1,3);

for i = 1:3
    e = zeros(3,1);
    e(i) = d;

    if cfg.meas.enable_lighttime
        fp = compute_rhoddot_lighttime(r_r + e, v_r, a_r, r_s, v_s, a_s, cfg.c);
        fm = compute_rhoddot_lighttime(r_r - e, v_r, a_r, r_s, v_s, a_s, cfg.c);
    else
        fp = compute_rhoddot(r_r + e, v_r, a_r, r_s, v_s, a_s);
        fm = compute_rhoddot(r_r - e, v_r, a_r, r_s, v_s, a_s);
    end

    Hpos(i) = (fp - fm) / (2*d);
end
end
