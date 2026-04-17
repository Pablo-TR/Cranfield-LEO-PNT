function rhoddot = compute_rhoddot(r_r, v_r, a_r, r_s, v_s, a_s)
dr  = r_s - r_r;
rho = norm(dr);
u   = dr / rho;

dv = v_s - v_r;
da = a_s - a_r;

rhodot  = dv.' * u;
rhoddot = da.' * u + (dv.'*dv - rhodot^2) / rho;
end
