function [pe, ve] = eci2ecef(pi_, vi_, utc)
%ECI2ECEF Convert inertial position/velocity to ECEF using GMST rotation.

we = 7.2921159e-5;
jd = juliandate_local(utc);
T  = (jd - 2451545) / 36525;
g  = 280.46061837 + 360.98564736629*(jd - 2451545) + ...
     0.000387933*T^2 - T^3/38710000;
th = deg2rad(mod(g, 360));

ct = cos(th);
st = sin(th);

R_ = [ct, st, 0; ...
     -st, ct, 0; ...
       0,  0, 1];

pe = R_ * pi_;

dR = we * [-st,  ct, 0; ...
           -ct, -st, 0; ...
             0,   0, 0];

ve = R_ * vi_ + dR * pi_;
end

function jd = juliandate_local(utc)
%JULIANDATE_LOCAL Toolbox-free Julian date for datetime inputs.

utc.TimeZone = 'UTC';
y = year(utc);
m = month(utc);
d = day(utc);

h = hour(utc);
mn = minute(utc);
s = second(utc);

idx = m <= 2;
y(idx) = y(idx) - 1;
m(idx) = m(idx) + 12;

A = floor(y/100);
B = 2 - A + floor(A/4);

day_fraction = (h + mn/60 + s/3600) / 24;
jd = floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + ...
     d + day_fraction + B - 1524.5;
end
