%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        mainMonteCarlo.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Analysis
% Author: Pablo
% Date: 12th April 2026
% Last Update: 12th April 2026
% Version: 1.0.0
%
% Description: Simple Monte Carlo study where each scenario is a random
% draw of all measurement-impairment parameters inside user-defined ranges.
% Several runs are executed per scenario, and the script stores flat
% vectors with all results together while printing global mean metrics.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

% Fix the random generator for repeatable studies.
rng(7);

% Get the current script folder.
scriptDir = fileparts(mfilename('fullpath'));
% Get the Tracking folder.
trackingDir = fileparts(scriptDir);
% Get the Algorithm folder.
algorithmDir = fileparts(trackingDir);
% Get the repository root folder.
repoRoot = fileparts(algorithmDir);
% Build the paths to the current tracking helpers.
kfHelperDir = fullfile(trackingDir, 'KFs and helpers');
simulatorDir = fullfile(trackingDir, 'Simulator');
% Add the tracking helpers, filters, and simulators to the MATLAB path.
addpath(kfHelperDir);
addpath(simulatorDir);

%% ========================================================================
%  STUDY CONFIGURATION
%  ========================================================================
% Set the number of random parameter combinations.
numScenarios = 50;
% Set the number of Monte Carlo runs per scenario.
runsPerScenario = 20;
% Set the Doppler sampling interval.
stepTime = 1e-3;
% Set the number of epochs per run.
numEpochs = 1000;
% Set the Doppler bin width used by the measurement generator.
binWidth = 500;
% Set the initial Doppler uncertainty for the filters.
sigmaFd = 450;
% Set the initial Doppler-rate uncertainty for the filters.
sigmaFddot = 5e4;
% Set the Doppler measurement standard deviation used by the filters.
measurementVar = 1000;
% Set the process-noise scale shared by the filters.
processNoiseScale = 800;
% Set the Doppler acceptance band used to declare convergence.
convergenceThresholdHz = 100;
% Require the Doppler error to stay inside the acceptance band for
% 100 consecutive samples before declaring convergence.
sustainSamples = 100;

% Enable the same Kalman filters used in the first Monte Carlo script.
runFilters = struct( ...
    'noVel', true, ...
    'velTle', true, ...
    'velTle2', true, ...
    'vel', true, ...
    'acc', true);

% Define short labels for the enabled filters.
filterLabels = struct( ...
    'noVel', '2nd order KF', ...
    'velTle', '2nd order KF, TLE velocity', ...
    'velTle2', '2nd order KF, error state', ...
    'vel', '2nd order KF, numerical velocity', ...
    'acc', '3rd order KF');

% Define the random ranges for each perturbation parameter.
perturbationRanges = struct( ...
    'neighbourJumpProbability', [0, 0.40], ...
    'outlierProbability', [0, 0.20], ...
    'outlierMeanHz', [-1000, 1000], ...
    'outlierStdHz', [100, 200], ...
    'gapProbability', [0.00], ...
    'oscillatorOffsetHz', [-50, 50], ...
    'oscillatorDriftRateHzPerSec', [-5, 5]);

% Define the TLE-based truth-generation configuration.
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

% Compute the total number of runs across all scenarios.
totalRuns = numScenarios * runsPerScenario;
% Create a flat counter for the global run index.
globalRunIndex = 0;

% Preallocate flat vectors with all sampled perturbation values.
flatParameters = struct( ...
    'neighbourJumpProbability', NaN(totalRuns, 1), ...
    'outlierProbability', NaN(totalRuns, 1), ...
    'outlierMeanHz', NaN(totalRuns, 1), ...
    'outlierStdHz', NaN(totalRuns, 1), ...
    'gapProbability', NaN(totalRuns, 1), ...
    'oscillatorOffsetHz', NaN(totalRuns, 1), ...
    'oscillatorDriftRateHzPerSec', NaN(totalRuns, 1));

% Preallocate the flat vector with the measurement RMSE of every run.
flatResults = struct('measurementRmse', NaN(totalRuns, 1));

% Preallocate the flat vectors with the filter metrics of every run.
filterNames = fieldnames(runFilters);
for filterIdx = 1:numel(filterNames)
    % Read the current filter name.
    filterName = filterNames{filterIdx};
    % Skip filters that are disabled in the configuration.
    if ~runFilters.(filterName)
        continue
    end

    % Create one flat vector per metric for the current filter.
    flatResults.(filterName) = struct( ...
        'rmse', NaN(totalRuns, 1), ...
        'postConvergenceRmse', NaN(totalRuns, 1), ...
        'convergenceTime', NaN(totalRuns, 1), ...
        'fdDotRmse', NaN(totalRuns, 1), ...
        'fdDotPostConvergenceRmse', NaN(totalRuns, 1));
end

%% ========================================================================
%  MONTE CARLO
%  ========================================================================
for scenarioIdx = 1:numScenarios
    % Draw one random perturbation configuration for this scenario.
    scenarioParameters = struct( ...
        'neighbourJumpProbability', sampleRange( ...
            perturbationRanges.neighbourJumpProbability), ...
        'outlierProbability', sampleRange( ...
            perturbationRanges.outlierProbability), ...
        'outlierMeanHz', sampleRange( ...
            perturbationRanges.outlierMeanHz), ...
        'outlierStdHz', sampleRange( ...
            perturbationRanges.outlierStdHz), ...
        'gapProbability', sampleRange( ...
            perturbationRanges.gapProbability), ...
        'oscillatorOffsetHz', sampleRange( ...
            perturbationRanges.oscillatorOffsetHz), ...
        'oscillatorDriftRateHzPerSec', sampleRange( ...
            perturbationRanges.oscillatorDriftRateHzPerSec));

    for runIdx = 1:runsPerScenario
        % Move to the next global run slot.
        globalRunIndex = globalRunIndex + 1;
        % Build a deterministic seed for this scenario/run pair.
        randomSeed = 100000 * scenarioIdx + runIdx;

        % Generate one Doppler measurement sequence for the current run.
        measurementData = generateDopplerMeasurements( ...
            'numEpochs', numEpochs, ...
            'stepTime', stepTime, ...
            'truthModel', 'tle', ...
            'truthParams', truthParams, ...
            'gaussianNoiseStdHz', measurementVar, ...
            'dopplerBinWidthHz', binWidth, ...
            'neighbourJumpProbability', ...
                scenarioParameters.neighbourJumpProbability, ...
            'outlierProbability', ...
                scenarioParameters.outlierProbability, ...
            'outlierMeanHz', scenarioParameters.outlierMeanHz, ...
            'outlierStdHz', scenarioParameters.outlierStdHz, ...
            'gapProbability', scenarioParameters.gapProbability, ...
            'oscillatorOffsetHz', scenarioParameters.oscillatorOffsetHz, ...
            'oscillatorDriftModel', 'linear', ...
            'oscillatorDriftParams', struct('driftRateHzPerSec', ...
                scenarioParameters.oscillatorDriftRateHzPerSec), ...
            'randomSeed', randomSeed);

        % Extract the truth time vector.
        time = measurementData.time(:);
        % Extract the true Doppler sequence.
        fdTrue = measurementData.fdTrue(:);
        % Extract the true Doppler-rate sequence.
        fdDotTrue = measurementData.fdDotTrue(:);
        % Prepare the measurements so filters can handle gaps safely.
        fdMeasForFilter = prepareMeasurementsForFilter( ...
            measurementData.fdMeasured);

        % Store the sampled perturbation values in flat vectors.
        flatParameters.neighbourJumpProbability(globalRunIndex) = ...
            scenarioParameters.neighbourJumpProbability;
        flatParameters.outlierProbability(globalRunIndex) = ...
            scenarioParameters.outlierProbability;
        flatParameters.outlierMeanHz(globalRunIndex) = ...
            scenarioParameters.outlierMeanHz;
        flatParameters.outlierStdHz(globalRunIndex) = ...
            scenarioParameters.outlierStdHz;
        flatParameters.gapProbability(globalRunIndex) = ...
            scenarioParameters.gapProbability;
        flatParameters.oscillatorOffsetHz(globalRunIndex) = ...
            scenarioParameters.oscillatorOffsetHz;
        flatParameters.oscillatorDriftRateHzPerSec(globalRunIndex) = ...
            scenarioParameters.oscillatorDriftRateHzPerSec;

        % Store the raw measurement RMSE for this run.
        flatResults.measurementRmse(globalRunIndex) = trackingRmse( ...
            fdTrue, measurementData.fdMeasured);

        for filterIdx = 1:numel(filterNames)
            % Read the current filter name.
            filterName = filterNames{filterIdx};
            % Skip filters that are disabled in the configuration.
            if ~runFilters.(filterName)
                continue
            end

            % Run the current Kalman filter.
            filterOutput = runFilter(filterName, fdMeasForFilter, ...
                fdTrue, fdDotTrue, stepTime, sigmaFd, sigmaFddot, ...
                processNoiseScale, measurementVar);

            % Compute the convergence time of the filter output.
            convergenceMetric = trackingConvergenceTime( ...
                fdTrue, filterOutput.fdEst, time, ...
                convergenceThresholdHz, sustainSamples);
            % Compute the post-convergence RMSE of the filter output.
            postConvergenceMetric = trackingPostConvergenceRmse( ...
                fdTrue, filterOutput.fdEst, time, ...
                convergenceThresholdHz, sustainSamples);

            % Store the total RMSE of the filter output.
            flatResults.(filterName).rmse(globalRunIndex) = ...
                trackingRmse(fdTrue, filterOutput.fdEst);
            % Store the post-convergence RMSE of the filter output.
            flatResults.(filterName).postConvergenceRmse(globalRunIndex) = ...
                postConvergenceMetric.rmseValue;
            % Store the convergence time of the filter output.
            flatResults.(filterName).convergenceTime(globalRunIndex) = ...
                convergenceMetric.convergenceTime;
            % Store the Doppler-rate RMSE of the filter output.
            flatResults.(filterName).fdDotRmse(globalRunIndex) = ...
                trackingRmse(fdDotTrue, filterOutput.fdDotEst);
            % Store the post-convergence Doppler-rate RMSE using the
            % Doppler convergence index as the only convergence criterion.
            flatResults.(filterName).fdDotPostConvergenceRmse( ...
                globalRunIndex) = postConvergenceRmseFromIndex( ...
                fdDotTrue, filterOutput.fdDotEst, ...
                convergenceMetric.convergenceIndex);
        end
    end
end

%% ========================================================================
%  GLOBAL SUMMARY
%  ========================================================================
% Compute the global mean measurement RMSE over all scenarios and runs.
summaryResults.measurementMeanRmse = meanFinite( ...
    flatResults.measurementRmse);

fprintf('\nGlobal Monte Carlo summary across all scenarios and runs:\n');
fprintf('Measurements mean RMSE: %.3f Hz\n', ...
    summaryResults.measurementMeanRmse);

for filterIdx = 1:numel(filterNames)
    % Read the current filter name.
    filterName = filterNames{filterIdx};
    % Skip filters that are disabled in the configuration.
    if ~runFilters.(filterName)
        continue
    end

    % Compute the global mean RMSE for the current filter.
    summaryResults.(filterName).meanRmse = meanFinite( ...
        flatResults.(filterName).rmse);
    % Compute the global mean post-convergence RMSE for the current filter.
    summaryResults.(filterName).meanPostConvergenceRmse = ...
        meanFinite(flatResults.(filterName).postConvergenceRmse);
    % Compute the global mean convergence time for the current filter.
    summaryResults.(filterName).meanConvergenceTime = ...
        meanFinite(flatResults.(filterName).convergenceTime);
    % Compute the global mean Doppler-rate RMSE for the current filter.
    summaryResults.(filterName).meanFdDotRmse = meanFinite( ...
        flatResults.(filterName).fdDotRmse);
    % Compute the global mean post-convergence Doppler-rate RMSE.
    summaryResults.(filterName).meanFdDotPostConvergenceRmse = ...
        meanFinite(flatResults.(filterName).fdDotPostConvergenceRmse);

    % Print the global mean RMSE for the current filter.
    fprintf('%s mean RMSE: %.3f Hz\n', ...
        filterLabels.(filterName), ...
        summaryResults.(filterName).meanRmse);
    % Print the global mean post-convergence RMSE for the current filter.
    fprintf('%s mean post-convergence RMSE: %.3f Hz\n', ...
        filterLabels.(filterName), ...
        summaryResults.(filterName).meanPostConvergenceRmse);
    % Print the global mean convergence time for the current filter.
    fprintf('%s mean convergence time: %.6f s\n', ...
        filterLabels.(filterName), ...
        summaryResults.(filterName).meanConvergenceTime);
    % Print the global mean Doppler-rate RMSE for the current filter.
    fprintf('%s mean Doppler-rate RMSE: %.3f Hz/s\n', ...
        filterLabels.(filterName), ...
        summaryResults.(filterName).meanFdDotRmse);
    % Print the global mean post-convergence Doppler-rate RMSE.
    fprintf('%s mean post-convergence Doppler-rate RMSE: %.3f Hz/s\n', ...
        filterLabels.(filterName), ...
        summaryResults.(filterName).meanFdDotPostConvergenceRmse);
end

%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================
function sampledValue = sampleRange(valueRange)

    % Return the fixed value directly when only one value is provided.
    if numel(valueRange) == 1
        sampledValue = valueRange;
        return
    end

    % Draw a uniform random sample inside the requested range.
    sampledValue = valueRange(1) + ...
        rand * (valueRange(2) - valueRange(1));

end

function filterOutput = runFilter(filterName, fdMeasForFilter, fdTrue, ...
    fdDotTrue, stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
    measurementVar)

    % Route the measurement sequence to the requested Kalman filter.
    switch filterName
        case 'noVel'
            filterOutput = kfDopplerNoVel(fdMeasForFilter, ...
                stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
                measurementVar);

        case 'velTle'
            filterOutput = kfDopplerVel_tle(fdMeasForFilter, ...
                fdDotTrue, stepTime, sigmaFd, sigmaFddot, ...
                processNoiseScale, measurementVar);

        case 'velTle2'
            filterOutput = kfDopplerVel_tle2(fdMeasForFilter, ...
                fdTrue, stepTime, sigmaFd, sigmaFddot, ...
                processNoiseScale, fdDotTrue, measurementVar);

        case 'vel'
            filterOutput = kfDopplerVel(fdMeasForFilter, ...
                stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
                measurementVar);

        case 'acc'
            filterOutput = kfDopplerAcc(fdMeasForFilter, ...
                stepTime, sigmaFd, sigmaFddot, processNoiseScale, ...
                measurementVar);

        otherwise
            error('Unsupported filter name: %s', filterName);
    end

end

function fdPrepared = prepareMeasurementsForFilter(fdMeasured)

    % Convert the input measurements into a column vector.
    fdPrepared = fdMeasured(:);
    % Fill missing samples by linear interpolation.
    fdPrepared = fillmissing(fdPrepared, 'linear', 'EndValues', 'nearest');

    % Replace any remaining invalid samples with a finite fallback value.
    if any(~isfinite(fdPrepared))
        % Extract the finite samples that remain after interpolation.
        finiteSamples = fdPrepared(isfinite(fdPrepared));
        % Use zero if every sample is still invalid.
        if isempty(finiteSamples)
            fdPrepared(:) = 0;
        else
            % Use the first finite value as a fallback for invalid entries.
            fdPrepared(~isfinite(fdPrepared)) = finiteSamples(1);
        end
    end

end

function meanValue = meanFinite(values)

    % Keep only the finite samples of the input vector.
    finiteValues = values(isfinite(values));

    % Return NaN when no finite values are available.
    if isempty(finiteValues)
        meanValue = NaN;
        return
    end

    % Compute the arithmetic mean of the finite values.
    meanValue = mean(finiteValues);

end

function rmseValue = postConvergenceRmseFromIndex(reference, estimate, ...
    startIndex)

    % Return NaN when the convergence index is not available.
    if ~isfinite(startIndex)
        rmseValue = NaN;
        return
    end

    % Compute the RMSE from the Doppler convergence sample to the end.
    rmseValue = trackingRmse(reference(startIndex:end), ...
        estimate(startIndex:end));

end
%%

