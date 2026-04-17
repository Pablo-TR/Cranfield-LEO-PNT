function delay_m = tropo_delay_simple(elev_deg, ztd_m, el_floor, eps_deg)
el = max(elev_deg, el_floor);
delay_m = ztd_m ./ sind(el + eps_deg);
end
