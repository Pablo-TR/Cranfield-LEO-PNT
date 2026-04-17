function rhoddot = compute_rhoddot_lighttime(r_r, v_r, a_r, r_s, v_s, a_s, c)
dr0    = r_s - r_r;
tau    = norm(dr0)/c;

r_s_tx = r_s - v_s*tau + 0.5*a_s*tau^2;
v_s_tx = v_s - a_s*tau;

rhoddot = compute_rhoddot(r_r, v_r, a_r, r_s_tx, v_s_tx, a_s);
end
