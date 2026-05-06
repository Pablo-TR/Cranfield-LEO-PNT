%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         mainVisualiseMeasurementPerturbations.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Analysis
% Author: Pablo
% Date: 24th March 2026
% Last Update: 24th March 2026
% Version: 1.0.0
%
% Description: Simple script to visualise the effect of each perturbation
% available in the Doppler measurement generator. Only one perturbation is
% enabled at a time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Clear the workspace variables.
clear;
% Close all open figures from previous runs.
close all;

% Force LaTeX formatting for figure text.
set(groot, 'defaultTextInterpreter', 'latex');
% Force LaTeX formatting for axis tick labels.
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
% Force LaTeX formatting for legends.
set(groot, 'defaultLegendInterpreter', 'latex');
% Set the default font size used in axes.
set(groot, 'defaultAxesFontSize', 14);
% Set the default font size used in text objects.
set(groot, 'defaultTextFontSize', 14);
% Set the default font size used in legends.
set(groot, 'defaultLegendFontSize', 14);

% Resolve the folder that contains this script.
scriptDir = fileparts(mfilename('fullpath'));
% Resolve the parent Tracking folder.
trackingDir = fileparts(scriptDir);
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
% Set the number of epochs generated in each perturbation example.
numEpochs = 1000;
% Define a simple polynomial truth model kept as an alternative reference.
truthParams = struct( ...
    'fd0', 4000, ...
    'fdDot', -1500, ...
    'fdDdot', 5000);

% Define the TLE-driven truth-model configuration used in the figures.
truthParams2 = struct( ...
    'tleFile', 'C:\Users\padlo\Desktop\Cranfield-LEO-PNT\Signals\oneweb.tle', ...
    'rxLatDeg', 52.06869, ...
    'rxLonDeg', -0.60884, ...
    'rxAltM', 100, ...
    'carrierFreqHz', 2.9e9, ...
    'startFromTLE', 'first', ...
    'selectionMode', 'highestMaxElevation', ...
    'minElevDeg', 60, ...
    'minMaxElevDeg', 70, ...
    'coarseStep', 10, ...
    'useSGP4', true, ...
    'fd0',0, ...
    'fdDot', 5, ...
    'fdDdot', 0.2, ...
    'allowJ2Fallback', true);

% Build the baseline argument list used by all perturbation examples.
baseArgs = { ...
    'numEpochs', numEpochs, ...
    'stepTime', stepTime, ...
    'truthModel', 'tle', ...
    'truthParams', truthParams2, ...
    'gaussianNoiseStdHz', 0, ...
    'dopplerBinWidthHz', 0, ...
    'neighbourJumpProbability', 0, ...
    'outlierProbability', 0, ...
    'outlierMeanHz', 0, ...
    'outlierStdHz', 0, ...
    'gapProbability', 0, ...
    'oscillatorOffsetHz', 0, ...
    'oscillatorDriftModel', 'linear', ...
    'oscillatorDriftParams', struct('driftRateHzPerSec', 0), ...
    'randomSeed', 7};

% Generate the reference sequence without perturbations.
data = generateDopplerMeasurements(baseArgs{:});
% Plot the reference sequence without perturbations.
plotPerturbation(data, 'Reference without perturbation');

% Generate the sequence with Gaussian noise only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'gaussianNoiseStdHz', 50);
% Plot the Gaussian-noise example.
plotPerturbation(data, 'Gaussian noise');

% Generate the sequence with Doppler-bin quantisation only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'dopplerBinWidthHz', 25);
% Plot the bin-quantisation example.
plotPerturbation(data, 'Doppler bin quantisation');

% Generate the sequence with neighbour-bin jumps on top of bin quantisation.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'dopplerBinWidthHz', 25, ...
    'neighbourJumpProbability', 0.05);
% Plot the neighbour-bin jump example.
plotPerturbation(data, 'Neighbour-bin jumps');

% Generate the sequence with outliers only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'outlierProbability', 0.05, ...
    'outlierMeanHz', 0, ...
    'outlierStdHz', 200);
% Plot the outlier example.
plotPerturbation(data, 'Outliers');

% Generate the sequence with measurement gaps only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'gapProbability', 0.08);
% Plot the measurement-gap example.
plotPerturbation(data, 'Measurement gaps');

% Generate the sequence with oscillator offset only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'oscillatorOffsetHz', 600);
% Plot the oscillator-offset example.
plotPerturbation(data, 'Oscillator offset');

% Generate the sequence with slow oscillator drift only.
data = generateDopplerMeasurements( ...
    baseArgs{:}, ...
    'oscillatorDriftModel', 'linear', ...
    'oscillatorDriftParams', struct('driftRateHzPerSec', 1500));
% Plot the oscillator-drift example.
plotPerturbation(data, 'Oscillator drift');

function plotPerturbation(data, plotTitleText)

    % Create a new figure named after the perturbation being visualised.
    figure('Name', plotTitleText);
    % Plot the underlying Doppler truth.
    plot(data.time, data.fdTrue, 'k', 'LineWidth', 2.5, ...
        'DisplayName', 'Ground truth');
    % Keep the figure active so the measurements can be overlaid.
    hold on;
    % Plot the degraded measurements created by the perturbation.
    plot(data.time, data.fdMeasured, 'b.-', 'LineWidth', 2.0, ...
        'DisplayName', 'Measured');
    % Label the horizontal axis.
    xlabel('$t$ [s]');
    % Label the vertical axis.
    ylabel('$f_D$ [Hz]');
    % Add the title of the perturbation figure.
    title(plotTitleText);
    % Add the legend that identifies truth and measurements.
    legend('Location', 'best');
    % Enable the plot grid.
    grid on;

end

