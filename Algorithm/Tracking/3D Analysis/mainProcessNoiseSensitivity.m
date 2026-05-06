%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              mainProcessNoiseSensitivity.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Analysis
% Author: Pablo
% Date: 15th March 2026
% Last Update: 15th March 2026
% Version: 1.0.0
%
% Description: Minimal sensitivity study of the Doppler tracking Kalman
% filters with respect to the process noise scale. The script evaluates
% convergence time, total RMSE, and post-convergence RMSE for the Doppler
% frequency estimate and plots each metric separately as a function of the
% process noise.
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

% Resolve the folder that contains this script.
scriptDir = fileparts(mfilename('fullpath'));
% Resolve the parent Tracking folder.
trackingDir = fileparts(scriptDir);
% Resolve the folder that stores the Kalman filters and helper metrics.
kfHelperDir = fullfile(trackingDir, 'KFs and helpers');
% Resolve the folder that stores the Doppler simulators.
simulatorDir = fullfile(trackingDir, 'Simulator');
% Add the Kalman filters and helper metrics to the MATLAB path.
addpath(kfHelperDir);
% Add the simulator functions to the MATLAB path.
addpath(simulatorDir);

%% ========================================================================
%  STUDY CONFIGURATION
%  ========================================================================
% Set the measurement sampling interval in seconds.
stepTime = 1e-3;
% Set the number of simulated Doppler samples.
numEpochs = 1000;
% Set the Doppler search-bin width used by the measurement generator.
binWidth = 350;
% Set the initial Doppler uncertainty used by the Kalman filters.
sigmaFd = 450;
% Set the initial Doppler-rate uncertainty used by the Kalman filters.
sigmaFddot = 1e7;
% Set the Doppler measurement noise level in Hz.
measurementVar = 250;
% Set the Doppler acceptance band used to declare convergence.
convergenceThresholdHz = 100;
% Set the number of consecutive samples required for convergence.
sustainSamples = 10;
% Create the logarithmic sweep of process-noise values to be tested.
processNoiseValues = logspace(-3, 3);

% Select which Kalman filters are included in the sensitivity study.
runFilters = struct( ...
    'noVel', true, ...
    'velTle', true, ...
    'velTle2', true, ...
    'vel', true, ...
    'acc', true);

% Define the legend labels associated with each filter.
filterLabels = struct( ...
    'noVel', 'No velocity KF', ...
    'velTle', 'Velocity TLE KF', ...
    'velTle2', 'Error-state TLE KF', ...
    'vel', 'Velocity KF', ...
    'acc', 'Acceleration KF');
% Initialise the structure that stores the study results.
results = struct();
% Read the list of configured filter names.
filterNames = fieldnames(runFilters);

% Preallocate one result entry for each enabled filter.
for fieldIndex = 1:numel(filterNames)
    % Read the current filter name.
    filterName = filterNames{fieldIndex};
    % Skip filters that are disabled in the configuration.
    if runFilters.(filterName)
        % Create the vectors that store the metrics of the current filter.
        results.(filterName) = struct( ...
            'processNoise', processNoiseValues(:), ...
            'convergenceTime', NaN(numel(processNoiseValues), 1), ...
            'totalRmse', NaN(numel(processNoiseValues), 1), ...
            'postConvergenceRmse', NaN(numel(processNoiseValues), 1), ...
            'finalDopplerVariance', NaN(numel(processNoiseValues), 1));
    end
end


%% ========================================================================
%  DOPPLER MEASUREMENT GENERATION
%  ========================================================================
% Define the TLE-driven truth-model configuration used for the study.
truthParams = struct( ...
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
    'allowJ2Fallback', true);

% Generate one Doppler measurement sequence with the nominal settings.
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
% Extract the true Doppler-rate sequence.
fdDotTrue = measurementData.fdDotTrue(:);
% Fill measurement gaps so the filters can process the full sequence.
fdMeasForFilter = prepareMeasurementsForFilter( ...
    measurementData.fdMeasured);

%% ========================================================================
%  PROCESS NOISE SWEEP
%  ========================================================================

% Sweep every candidate process-noise value.
for processNoiseIndex = 1:numel(processNoiseValues)
    % Read the current process-noise value.
    processNoiseScale = processNoiseValues(processNoiseIndex);

    % Run the no-velocity filter when it is enabled.
    if runFilters.noVel
        output = kfDopplerNoVel(fdMeasForFilter, stepTime, ...
            sigmaFd, sigmaFddot, processNoiseScale, measurementVar);
        % Store the metrics produced by the no-velocity filter.
        results.noVel = storeMetrics(results.noVel, processNoiseIndex, ...
            fdTrue, output, time, convergenceThresholdHz, ...
            sustainSamples);
    end

    % Run the TLE-velocity filter when it is enabled.
    if runFilters.velTle
        output = kfDopplerVel_tle(fdMeasForFilter, fdDotTrue, ...
            stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
            measurementVar);
        % Store the metrics produced by the TLE-velocity filter.
        results.velTle = storeMetrics(results.velTle, ...
            processNoiseIndex, fdTrue, output, time, ...
            convergenceThresholdHz, sustainSamples);
    end

    % Run the TLE error-state filter when it is enabled.
    if runFilters.velTle2
        output = kfDopplerVel_tle2(fdMeasForFilter, fdTrue, ...
            stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
            fdDotTrue, measurementVar);
        % Store the metrics produced by the TLE error-state filter.
        results.velTle2 = storeMetrics(results.velTle2, ...
            processNoiseIndex, fdTrue, output, time, ...
            convergenceThresholdHz, sustainSamples);
    end

    % Run the numerical-velocity filter when it is enabled.
    if runFilters.vel
        output = kfDopplerVel(fdMeasForFilter, stepTime, sigmaFd, ...
            sigmaFddot, processNoiseScale, measurementVar);
        % Store the metrics produced by the numerical-velocity filter.
        results.vel = storeMetrics(results.vel, processNoiseIndex, ...
            fdTrue, output, time, convergenceThresholdHz, ...
            sustainSamples);
    end

    % Run the acceleration-state filter when it is enabled.
    if runFilters.acc
        output = kfDopplerAcc(fdMeasForFilter, stepTime, ...
            sigmaFd, sigmaFddot, processNoiseScale, measurementVar);
        % Store the metrics produced by the acceleration-state filter.
        results.acc = storeMetrics(results.acc, processNoiseIndex, ...
            fdTrue, output, time, convergenceThresholdHz, ...
            sustainSamples);
    end
end

%% ========================================================================
%  PLOTS
%  ========================================================================
plotMetrics(processNoiseValues, results, runFilters, filterLabels, ...
    'convergenceTime', 'Convergence time [s]', ...
    'Convergence time vs process noise');

plotMetrics(processNoiseValues, results, runFilters, filterLabels, ...
    'totalRmse', 'Total RMSE [Hz]', ...
    'Total Doppler RMSE vs process noise');

plotMetrics(processNoiseValues, results, runFilters, filterLabels, ...
    'postConvergenceRmse', 'Post-convergence RMSE [Hz]', ...
    'Post-convergence Doppler RMSE vs process noise');

plotMetrics(processNoiseValues, results, runFilters, filterLabels, ...
    'finalDopplerVariance', 'Final Doppler variance [Hz$^2$]', ...
    'Final Doppler variance vs process noise');

%% ========================================================================
%   FUNCTIONS
%  ========================================================================
function fdPrepared = prepareMeasurementsForFilter(fdMeasured)

    % Convert the measurement sequence into a column vector.
    fdPrepared = fdMeasured(:);
    % Fill interior and edge gaps by linear interpolation.
    fdPrepared = fillmissing(fdPrepared, 'linear', 'EndValues', 'nearest');

    % Check whether any invalid samples remain after interpolation.
    if any(~isfinite(fdPrepared))
        % Keep only the finite samples that are still available.
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

function resultEntry = storeMetrics(resultEntry, processNoiseIndex, ...
    fdTrue, output, time, threshold, sustainSamples)

    % Extract the Doppler estimate produced by the current filter.
    fdEst = output.fdEst;

    % Compute the convergence metric of the Doppler estimate.
    convergenceMetric = trackingConvergenceTime(fdTrue, fdEst, ...
        time, threshold, sustainSamples);
    % Compute the post-convergence RMSE of the Doppler estimate.
    postConvergenceMetric = trackingPostConvergenceRmse(fdTrue, ...
        fdEst, time, threshold, sustainSamples);

    % Store the convergence time at the current process-noise index.
    resultEntry.convergenceTime(processNoiseIndex) = ...
        convergenceMetric.convergenceTime;
    % Store the total Doppler RMSE at the current process-noise index.
    resultEntry.totalRmse(processNoiseIndex) = ...
        trackingRmse(fdTrue, fdEst);
    % Store the post-convergence RMSE at the current process-noise index.
    resultEntry.postConvergenceRmse(processNoiseIndex) = ...
        postConvergenceMetric.rmseValue;
    % Store the final Doppler variance of the filter state.
    resultEntry.finalDopplerVariance(processNoiseIndex) = ...
        output.covarianceHistory(1, 1, end);

end

function plotMetrics(processNoiseValues, results, runFilters, ...
    filterLabels, metricField, yLabelText, figureTitleText)

    % Create a new figure for the current metric.
    figure;
    % Keep the figure active while the filter curves are added.
    hold on;

    % Plot the no-velocity filter curve when it is enabled.
    if runFilters.noVel
        semilogx(processNoiseValues, results.noVel.(metricField), ...
            'b-o', 'LineWidth', 1.2, 'DisplayName', filterLabels.noVel);
    end

    % Plot the TLE-velocity filter curve when it is enabled.
    if runFilters.velTle
        semilogx(processNoiseValues, results.velTle.(metricField), ...
            'Color', [0.85 0.45 0.10], 'Marker', 's', ...
            'LineWidth', 1.2, 'DisplayName', filterLabels.velTle);
    end

    % Plot the TLE error-state filter curve when it is enabled.
    if runFilters.velTle2
        semilogx(processNoiseValues, results.velTle2.(metricField), ...
            'Color', [0.10 0.60 0.70], 'Marker', 'd', ...
            'LineWidth', 1.2, 'DisplayName', filterLabels.velTle2);
    end

    % Plot the numerical-velocity filter curve when it is enabled.
    if runFilters.vel
       semilogx(processNoiseValues, results.vel.(metricField), ...
            'm-^', 'LineWidth', 1.2, 'DisplayName', filterLabels.vel);
    end

    % Plot the acceleration-state filter curve when it is enabled.
    if runFilters.acc
        semilogx(processNoiseValues, results.acc.(metricField), ...
            'r-v', 'LineWidth', 1.2, 'DisplayName', filterLabels.acc);
    end

    % Label the horizontal axis.
    xlabel('Process noise scale');
    % Label the vertical axis.
    ylabel(yLabelText);
    % Add the figure title.
    title(figureTitleText);
    % Add the legend that identifies each filter curve.
    legend('Location', 'best');
    % Enable the plot grid.
    grid on;

end

