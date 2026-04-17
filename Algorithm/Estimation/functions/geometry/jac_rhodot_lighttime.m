function H6 = jac_rhodot_lighttime(r_r, v_r, r_s, v_s, c)
dr0    = r_s - r_r;
tau    = norm(dr0)/c;
r_s_tx = r_s - v_s*tau;
H6     = jac_rhodot(r_r, v_r, r_s_tx, v_s);
end
