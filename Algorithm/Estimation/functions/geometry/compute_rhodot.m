function rhodot = compute_rhodot(r_r, v_r, r_s, v_s)
dr = r_s - r_r;
u  = dr / norm(dr);
dv = v_s - v_r;
rhodot = dv.' * u;
end
