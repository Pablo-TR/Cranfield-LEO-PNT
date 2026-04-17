function el = compute_elev(spos, rpos, R_e2n)
drENU = R_e2n * (spos - rpos);
el    = rad2deg(atan2(drENU(3), hypot(drENU(1),drENU(2))));
end
