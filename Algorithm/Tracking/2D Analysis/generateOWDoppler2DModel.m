%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               generateDoppler2DModel.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 4th May 2026
% Last Update: 5th May 2026
% Version: 1.0.0
%
% Description: Ground-truth generator based on a circular orbit, circular
% Earth model. Satellite is passing overhead at a LEO altitude
% corresponding to OneWeb altitude. Earth is not rotating. Assumes receiver
% at Earth's surface and no elevation masking.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,fd, dfdt, d2fdt2] = generateOWDoppler2DModel(timestep,plotFigure, h, f0, fmix, lambda0)
% Replace the sentinel value with the default time step.
if timestep == -1
    timestep = 2.5;
end

% Build the simulation time axis.
t = 0:timestep:1300;
% Replace the sentinel value with the default orbital height.
if h == -1
    h  = 1200e3;        % Orbital height [m]
end

% Define the Earth radius.
Re = 6.3781e6;      % Earth radius [m]
% Define the Earth gravitational parameter.
mu = 3.986e14;      % Earth's mu [m^3/s^2]
% Define the speed of light.
c  = 3e8;           % Speed of light [m/s]
% Replace the sentinel value with the default carrier frequency.
if f0 == -1
    f0 = 11.325e9;         % Carrier [Hz]
end

% Replace the sentinel value with the default local-oscillator frequency.
if fmix == -1
    fmix = f0;          % LO [Hz]
end
% Compute the orbital radius.
r     = Re + h;
% Compute the orbital speed.
Vorb  = sqrt(mu/r);
% Compute the orbital angular rate.
omega = sqrt(mu/r^3);

% lambda = central angle wrt receiver zenith
% Replace the sentinel value with the default initial central angle.
if lambda0 == -1
    lambda0 = -0.6;                 % start before closest approach [rad]
end

% Propagate the central angle over time.
lambda  = lambda0 + omega*t;

% Slant range
% Compute the slant range from the receiver to the satellite.
rho = sqrt(r^2 + Re^2 - 2*r*Re*cos(lambda));

% Radial velocity: positive when approaching
% Compute the line-of-sight velocity.
Vrad = -Vorb .* Re .* sin(lambda) ./ rho;

% Received frequency and Doppler
% Convert the radial velocity into received carrier frequency.
frec = f0 .* (1 + Vrad./c);
% Convert the received frequency into Doppler shift.
fd   = frec - fmix;
% Convert the Doppler vector into a column vector.
fd = fd';

% Obtain maximum Doppler rate and maximum Doppler acceleration
% Compute the Doppler rate numerically.
dfdt = diff(fd)/timestep;
% Compute the Doppler acceleration numerically.
d2fdt2 = diff(dfdt)/timestep;

% Keep only the maximum absolute Doppler acceleration.
d2fdt2 = max(abs(d2fdt2));


% Plot the simulated Doppler curve when requested.
if plotFigure
    plot(t, fd.*1e-3, 'LineWidth', 1.5)
    grid on
    xlabel('Time [s]')
    ylabel('Doppler [kHz]')
    title('Doppler during 2D overhead pass')
end
