function [x,y,z] = geodetic2ecef_local(lat_deg, lon_deg, alt_m)
a   = 6378137.0;
e2  = 0.00669437999014;
lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);
N   = a / sqrt(1 - e2*sin(lat)^2);

x = (N+alt_m)*cos(lat)*cos(lon);
y = (N+alt_m)*cos(lat)*sin(lon);
z = (N*(1-e2)+alt_m)*sin(lat);
end
