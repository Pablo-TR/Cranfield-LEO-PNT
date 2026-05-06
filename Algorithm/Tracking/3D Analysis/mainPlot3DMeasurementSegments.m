%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              mainPlot3DMeasurementSegments.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Analysis
% Author: Pablo
% Date: 15th April 2026
% Last Update: 15th April 2026
% Version: 1.0.0
%
% Description: Simple script that uses the TLE-based 3D Doppler simulator
% to generate noisy measurements and plot the full ground-truth Doppler
% curve together with a few short measured segments distributed along it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Clear the workspace variables.
clear;
% Close all open figures from previous runs.
close all;

% Fix the random generator for repeatable plots.
rng(7);

% Resolve the folder that contains this script.
scriptDir = fileparts(mfilename('fullpath'));
% Resolve the parent Tracking folder.
trackingDir = fileparts(scriptDir);
% Resolve the Algorithm folder.
algorithmDir = fileparts(trackingDir);
% Resolve the repository root folder.
repoRoot = fileparts(algorithmDir);
% Resolve the folder that stores the Kalman filters and helper metrics.
kfHelperDir = fullfile(trackingDir, 'KFs and helpers');
% Resolve the folder that stores the simulator functions.
simulatorDir = fullfile(trackingDir, 'Simulator');
% Add the Kalman filters and helper metrics to the MATLAB path.
addpath(kfHelperDir);
% Add the simulator functions to the MATLAB path.
addpath(simulatorDir);

% Set the sampling interval of the generated measurements in seconds.
stepTime = 1;
% Set the number of epochs generated in the plot.
numEpochs = 1000;
% Set the standard deviation of the measurement noise in Hz.
measurementNoiseStdHz = 10000;
% Store the nominal segment length used as a visual reference.
segmentLength = 600;
% Store the number of short measurement segments overlaid on the truth.
numberOfSegments = 4;

% Define the TLE-driven truth-model configuration used in the figure.
truthParams = struct( ...
    'tleFile', fullfile(repoRoot, 'Signals', 'oneweb.tle'), ...
    'rxLatDeg', 52.06869, ...
    'rxLonDeg', -0.60884, ...
    'rxAltM', 100, ...
    'carrierFreqHz', 11.325e9, ...
    'startFromTLE', 'first', ...
    'selectionMode', 'highestMaxElevation', ...
    'minElevDeg', 60, ...
    'minMaxElevDeg', 70, ...
    'coarseStep', 10, ...
    'useSGP4', true, ...
    'allowJ2Fallback', true);

% Generate one TLE-based Doppler sequence with additive measurement noise.
measurementData = generateDopplerMeasurements( ...
    'numEpochs', numEpochs, ...
    'stepTime', stepTime, ...
    'truthModel', 'tle', ...
    'truthParams', truthParams, ...
    'gaussianNoiseStdHz', measurementNoiseStdHz, ...
    'randomSeed', 7);

% Extract the simulation time vector as a column vector.
time = measurementData.time(:);
% Convert the Doppler truth from Hz to kHz for plotting.
fdTruekHz = measurementData.fdTrue(:) / 1e3;
% Convert the noisy Doppler measurements from Hz to kHz for plotting.
fdMeasuredkHz = measurementData.fdMeasured(:) / 1e3;

%%
% Force LaTeX formatting for figure text.
set(groot, 'defaultTextInterpreter', 'latex');
% Force LaTeX formatting for axis tick labels.
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% Force LaTeX formatting for legends.
set(groot, 'defaultLegendInterpreter', 'latex');
% Set the default font size used in axes.
set(groot, 'defaultAxesFontSize', 20);
% Set the default font size used in text objects.
set(groot, 'defaultTextFontSize', 14);
% Set the default font size used in legends.
set(groot, 'defaultLegendFontSize', 14);
% Create the figure that shows the global S-curve and the selected segments.
figure('Name', '3D Doppler curve with measured segments');
% Plot the full ground-truth Doppler curve.
plot(time, fdTruekHz, '--', 'LineWidth', 1.5, ...
    'DisplayName', 'Ground truth');
% Keep the figure active so the measurement segments can be overlaid.
hold on;

% Plot the first short measurement segment.
plot(time(200:300), fdMeasuredkHz(200:300), 'LineWidth', 1.5, 'Color','red');
% Plot the second short measurement segment.
plot(time(450:550), fdMeasuredkHz(450:550), 'LineWidth', 1.5, 'Color','red');
% Plot the third short measurement segment.
plot(time(600:700), fdMeasuredkHz(600:700), 'LineWidth', 1.5, 'Color','red');
% Plot the fourth short measurement segment.
plot(time(350:375), fdMeasuredkHz(350:375), 'LineWidth', 1.5, 'Color','red');

% Enable the plot grid.
grid on;
% Label the horizontal axis.
xlabel('Time [s]', 'FontSize',20);
% Label the vertical axis.
ylabel('Doppler [kHz]', 'FontSize',20);

