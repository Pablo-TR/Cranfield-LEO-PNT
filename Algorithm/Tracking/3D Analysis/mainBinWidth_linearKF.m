%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                mainBinWidth_linearKF.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Analysis
% Author: Pablo
% Date: 16th April 2026
% Last Update: 16th April 2026
% Version: 1.0.0
%
% Description: Minimal sensitivity study of the Doppler tracking linearKF
% filters with respect to the Doppler-bin width. The script evaluates
% convergence time, total RMSE, and post-convergence RMSE for the Doppler
% frequency estimate.
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

% Fix the random generator for repeatable studies.
rng(7);

%% ========================================================================
%  STUDY CONFIGURATION
%  ========================================================================
% Set the Doppler sampling interval in seconds.
stepTime = 1;
% Set the number of Doppler samples in each simulation.
numEpochs = 1000;
% Set the measurement noise level in Hz.
measurementVar = 500;
% Set the convergence threshold in Hz.
convergenceThresholdHz = 100;
% Set the number of consecutive samples required for convergence.
sustainSamples = 5;
% Set the process-noise scale shared by the linear filters.
processNoiseScale = 10;
% Define the Doppler-bin widths explored by the study.
binWidthValues = 0:100:1000;

% Enable the linear Kalman-filter variants that will be compared.
runFilters = struct( ...
    'linear2State', true, ...
    'linear3State', true);

% Define the legend labels for the linear Kalman-filter variants.
filterLabels = struct( ...
    'linear2State', 'linearKF 2-state', ...
    'linear3State', 'linearKF 3-state');

% Initialise the results container.
results = struct();
% Read the list of configured filter names.
filterNames = fieldnames(runFilters);

% Preallocate one result entry per enabled filter.
for fieldIndex = 1:numel(filterNames)
    % Read the current filter name.
    filterName = filterNames{fieldIndex};
    % Create the result vectors only for enabled filters.
    if runFilters.(filterName)
        results.(filterName) = struct( ...
            'binWidth', binWidthValues(:), ...
            'convergenceTime', NaN(numel(binWidthValues), 1), ...
            'totalRmse', NaN(numel(binWidthValues), 1), ...
            'postConvergenceRmse', NaN(numel(binWidthValues), 1));
    end
end

%% ========================================================================
%  DOPPLER MEASUREMENT GENERATION
%  ========================================================================
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

% Define the TLE-driven truth-model configuration used in the study.
truthParams = struct( ...
    'tleFile', fullfile(repoRoot, 'Signals', 'oneweb.tle'), ...
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
    'allowJ2Fallback', true);

% Sweep every candidate Doppler-bin width.
for idx = 1:numel(binWidthValues)
    % Read the current Doppler-bin width.
    binWidth = binWidthValues(idx);

    % Generate one measurement sequence using the current bin width.
    measurementData = generateDopplerMeasurements( ...
        'numEpochs', numEpochs, ...
        'stepTime', stepTime, ...
        'truthModel', 'tle', ...
        'truthParams', truthParams, ...
        'gaussianNoiseStdHz', measurementVar, ...
        'dopplerBinWidthHz', binWidth, ...
        'neighbourJumpProbability', 0, ...
        'outlierProbability', 0, ...
        'gapProbability', 0, ...
        'oscillatorOffsetHz', 0, ...
        'oscillatorDriftModel', 'linear', ...
        'oscillatorDriftParams', struct('driftRateHzPerSec', 0), ...
        'randomSeed', 7);

    % Extract the simulation time vector as a column vector.
    time = measurementData.time(:);
    % Extract the true Doppler sequence.
    fdTrue = measurementData.fdTrue(:);
    % Extract the true Doppler acceleration sequence.
    fdDdotTrue = measurementData.fdDdotTrue(:);
    % Fill measurement gaps so the filters can process the full sequence.
    fdMeasForFilter = prepareMeasurementsForFilter( ...
        measurementData.fdMeasured);

    % Compute the maximum absolute Doppler acceleration used by linearKF.
    d2fdt2Max = max(abs(fdDdotTrue));
    % Replace invalid or zero acceleration bounds with a minimal fallback.
    if ~isfinite(d2fdt2Max) || d2fdt2Max <= 0
        d2fdt2Max = eps;
    end

    % Run the 2-state linear filter when it is enabled.
    if runFilters.linear2State
        output = linearKF(fdMeasForFilter, stepTime, processNoiseScale, ...
            measurementVar, d2fdt2Max, 0);
        % Store the metrics produced by the 2-state linear filter.
        results.linear2State = storeMetrics(results.linear2State, ...
            idx, fdTrue, output.fdEst, time, convergenceThresholdHz, ...
            sustainSamples);
    end

    % Run the 3-state linear filter when it is enabled.
    if runFilters.linear3State
        output = linearKF(fdMeasForFilter, stepTime, processNoiseScale, ...
            measurementVar, d2fdt2Max, 1);
        % Store the metrics produced by the 3-state linear filter.
        results.linear3State = storeMetrics(results.linear3State, ...
            idx, fdTrue, output.fdEst, time, convergenceThresholdHz, ...
            sustainSamples);
    end
end

%% ========================================================================
%  PLOTS
%  ========================================================================
plotMetrics(binWidthValues, results, runFilters, filterLabels, ...
    'convergenceTime', 'Convergence time [s]', ...
    'Convergence time vs bin width');

plotMetrics(binWidthValues, results, runFilters, filterLabels, ...
    'totalRmse', 'Total RMSE [Hz]', ...
    'Total Doppler RMSE vs bin width');

plotMetrics(binWidthValues, results, runFilters, filterLabels, ...
    'postConvergenceRmse', 'Post-convergence RMSE [Hz]', ...
    'Post-convergence Doppler RMSE vs bin width');

%% ========================================================================
%   LOCAL FUNCTIONS
%  ========================================================================
function fdPrepared = prepareMeasurementsForFilter(fdMeasured)

    % Convert the measurement sequence into a column vector.
    fdPrepared = fdMeasured(:);
    % Fill interior and edge gaps by linear interpolation.
    fdPrepared = fillmissing(fdPrepared, 'linear', 'EndValues', 'nearest');

    % Check whether any invalid samples remain after interpolation.
    if any(~isfinite(fdPrepared))
        % Keep only the finite samples that remain available.
        finiteSamples = fdPrepared(isfinite(fdPrepared));
        % Replace the whole sequence with zeros when nothing finite remains.
        if isempty(finiteSamples)
            fdPrepared(:) = 0;
        else
            % Replace the remaining invalid samples with the first finite value.
            fdPrepared(~isfinite(fdPrepared)) = finiteSamples(1);
        end
    end

end

function resultEntry = storeMetrics(resultEntry, idx, ...
    fdTrue, fdEst, time, threshold, sustainSamples)

    % Compute the convergence metric of the Doppler estimate.
    convergenceMetric = trackingConvergenceTime(fdTrue, fdEst, ...
        time, threshold, sustainSamples);
    % Compute the post-convergence RMSE of the Doppler estimate.
    postConvergenceMetric = trackingPostConvergenceRmse(fdTrue, ...
        fdEst, time, threshold, sustainSamples);

    % Store the convergence time at the current bin-width index.
    resultEntry.convergenceTime(idx) = ...
        convergenceMetric.convergenceTime;
    % Store the total Doppler RMSE at the current bin-width index.
    resultEntry.totalRmse(idx) = ...
        trackingRmse(fdTrue, fdEst);
    % Store the post-convergence RMSE at the current bin-width index.
    resultEntry.postConvergenceRmse(idx) = ...
        postConvergenceMetric.rmseValue;

end

function plotMetrics(values, results, runFilters, filterLabels, ...
    metricField, yLabelText, figureTitleText)

    % Create a new figure for the current metric.
    figure;
    % Keep the figure active while the filter curves are added.
    hold on;

    % Plot the 2-state linear filter curve when it is enabled.
    if runFilters.linear2State
        plot(values, results.linear2State.(metricField), ...
            'b-o', 'LineWidth', 1.2, ...
            'DisplayName', filterLabels.linear2State);
    end

    % Plot the 3-state linear filter curve when it is enabled.
    if runFilters.linear3State
        plot(values, results.linear3State.(metricField), ...
            'r-s', 'LineWidth', 1.2, ...
            'DisplayName', filterLabels.linear3State);
    end

    % Label the horizontal axis.
    xlabel('Bin Width [Hz]');
    % Label the vertical axis.
    ylabel(yLabelText);
    % Add the figure title.
    title(figureTitleText);
    % Add the legend that identifies each filter curve.
    legend('Location', 'best');
    % Enable the plot grid.
    grid on;

end
