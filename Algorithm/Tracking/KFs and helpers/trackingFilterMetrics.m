%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  trackingFilterMetrics.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.1.0
%
% Description: Summarises the main tracking metrics shown by the Kalman
% filter scripts, including RMSE, error extrema, convergence, sensitivity
% to initial error, and gap/outlier robustness.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metrics = trackingFilterMetrics(fdTrue, fdMeas, ...
    filterOutput, varargin)

    % Convert the true Doppler sequence into a column vector.
    fdTrue = fdTrue(:);
    % Convert the measured Doppler sequence into a column vector.
    fdMeas = fdMeas(:);

    % Extract the estimated Doppler sequence from the filter output.
    fdEst = filterOutput.fdEst(:);

    % Parse the optional truth vectors and configuration options.
    [fdDotTrue, fdDdotTrue, options] = parseOptionalInputs(varargin);
    % Resolve the time vector used by the convergence-based metrics.
    timeVector = resolveTimeVector(numel(fdTrue), filterOutput, options);
    % Resolve the Doppler convergence threshold used by the convergence-based metrics.
    threshold = resolveConvergenceThreshold(fdTrue, fdMeas, options);

    % Initialise the structure that will collect all tracking metrics.
    metrics = struct();
    % Compute the RMSE of the raw Doppler measurements.
    metrics.rmseMeasurementFd = trackingRmse(fdTrue, fdMeas);
    % Compute the RMSE of the estimated Doppler.
    metrics.rmseEstimatedFd = trackingRmse(fdTrue, fdEst);
    % Compute the maximum absolute error of the raw Doppler measurements.
    metrics.maxErrorMeasurementFd = trackingMaxError(fdTrue, fdMeas);
    % Compute the maximum absolute error of the estimated Doppler.
    metrics.maxErrorEstimatedFd = trackingMaxError(fdTrue, fdEst);
    % Compute the convergence metric of the raw Doppler measurements.
    metrics.convergenceMeasurementFd = trackingConvergenceTime( ...
        fdTrue, fdMeas, timeVector, threshold, ...
        options.convergenceSustainSamples);
    % Compute the convergence metric of the estimated Doppler.
    metrics.convergenceEstimatedFd = trackingConvergenceTime( ...
        fdTrue, fdEst, timeVector, threshold, ...
        options.convergenceSustainSamples);
    % Compute the post-convergence RMSE of the raw Doppler measurements.
    metrics.postConvergenceRmseMeasurementFd = ...
        trackingPostConvergenceRmse(fdTrue, fdMeas, timeVector, ...
        threshold, options.convergenceSustainSamples);
    % Compute the post-convergence RMSE of the estimated Doppler.
    metrics.postConvergenceRmseEstimatedFd = ...
        trackingPostConvergenceRmse(fdTrue, fdEst, timeVector, ...
        threshold, options.convergenceSustainSamples);
    % Compute the sensitivity of the estimated Doppler to the initial error.
    metrics.initialErrorSensitivityFd = ...
        trackingInitialErrorSensitivity(fdTrue, fdEst);

    % Compute the Doppler-rate RMSE when the truth and estimate are available.
    if ~isempty(fdDotTrue) && isfield(filterOutput, 'fdDotEst')
        metrics.rmseEstimatedFdDot = ...
            trackingRmse(fdDotTrue(:), filterOutput.fdDotEst(:));
    end

    % Compute the Doppler-acceleration metrics when the truth and estimate are available.
    if ~isempty(fdDdotTrue) && isfield(filterOutput, 'fdDdotEst')
        % Extract the estimated Doppler-acceleration sequence.
        fdDdotEst = filterOutput.fdDdotEst(:);
        % Compute the Doppler-acceleration RMSE.
        metrics.rmseEstimatedFdDdot = ...
            trackingRmse(fdDdotTrue(:), fdDdotEst);
        % Resolve the Doppler-acceleration convergence threshold.
        thresholdFdDdot = resolveAuxiliaryConvergenceThreshold( ...
            fdDdotTrue(:), fdDdotEst, ...
            options.convergenceThresholdFdDdotHzPerSec2);
        % Compute the Doppler-acceleration convergence metric.
        metrics.convergenceEstimatedFdDdot = ...
            trackingConvergenceTime(fdDdotTrue(:), fdDdotEst, ...
            timeVector, thresholdFdDdot, options.convergenceSustainSamples);
        % Compute the post-convergence RMSE of the Doppler-acceleration estimate.
        metrics.postConvergenceRmseEstimatedFdDdot = ...
            trackingPostConvergenceRmse(fdDdotTrue(:), fdDdotEst, ...
            timeVector, thresholdFdDdot, ...
            options.convergenceSustainSamples);
    end

    % Compute gap and outlier robustness only when masks are provided.
    if ~isempty(options.gapMask) || ~isempty(options.outlierMask)
        % Initialise the gap mask with a logical vector of false values.
        gapMask = false(size(fdTrue));
        % Initialise the outlier mask with a logical vector of false values.
        outlierMask = false(size(fdTrue));

        % Replace the default gap mask when the user provides one.
        if ~isempty(options.gapMask)
            gapMask = logical(options.gapMask(:));
        end
        % Replace the default outlier mask when the user provides one.
        if ~isempty(options.outlierMask)
            outlierMask = logical(options.outlierMask(:));
        end

        % Check that both masks have the same number of samples as the truth.
        if numel(gapMask) ~= numel(fdTrue) || ...
                numel(outlierMask) ~= numel(fdTrue)
            error(['gapMask and outlierMask must match the number of ' ...
                'samples in fdTrue.']);
        end

        % Compute the robustness of the raw Doppler measurements to gaps and outliers.
        metrics.gapOutlierRobustnessMeasurementFd = ...
            trackingGapOutlierRobustness(fdTrue, fdMeas, ...
            gapMask, outlierMask, options.recoveryWindowSamples);
        % Compute the robustness of the estimated Doppler to gaps and outliers.
        metrics.gapOutlierRobustnessEstimatedFd = ...
            trackingGapOutlierRobustness(fdTrue, fdEst, ...
            gapMask, outlierMask, options.recoveryWindowSamples);
    end

end

function [fdDotTrue, fdDdotTrue, options] = parseOptionalInputs(inputArgs)

    % Initialise the optional Doppler-rate truth as empty.
    fdDotTrue = [];
    % Initialise the optional Doppler-acceleration truth as empty.
    fdDdotTrue = [];
    % Build the default options structure used by the metric summary.
    options = struct( ...
        'timeVector', [], ...
        'stepTime', [], ...
        'convergenceThresholdHz', [], ...
        'convergenceThresholdFdDdotHzPerSec2', [], ...
        'convergenceSustainSamples', 5, ...
        'gapMask', [], ...
        'outlierMask', [], ...
        'recoveryWindowSamples', 5);

    % Count how many numeric optional inputs have been supplied.
    numericArgCount = 0;
    % Inspect every optional argument in order.
    for argIdx = 1:numel(inputArgs)
        % Read the current optional argument.
        currentArg = inputArgs{argIdx};
        % Merge fields when the current optional argument is a structure.
        if isstruct(currentArg)
            % Read the names of the fields stored in the structure.
            optionNames = fieldnames(currentArg);
            % Copy every supplied option field into the default structure.
            for nameIdx = 1:numel(optionNames)
                options.(optionNames{nameIdx}) = ...
                    currentArg.(optionNames{nameIdx});
            end
        else
            % Count the numeric optional arguments.
            numericArgCount = numericArgCount + 1;
            % Route the numeric optional arguments to the expected outputs.
            switch numericArgCount
                case 1
                    % Store the Doppler-rate truth vector.
                    fdDotTrue = currentArg(:);
                case 2
                    % Store the Doppler-acceleration truth vector.
                    fdDdotTrue = currentArg(:);
                otherwise
                    error('Too many numeric optional inputs were provided.');
            end
        end
    end

end

function threshold = resolveAuxiliaryConvergenceThreshold(reference, ...
    estimate, userThreshold)

    % Use the user-supplied threshold when it is provided explicitly.
    if ~isempty(userThreshold)
        threshold = userThreshold;
        return
    end

    % Compute the RMSE of the auxiliary estimate.
    estimateRmse = trackingRmse(reference(:), estimate(:));
    % Convert the RMSE into a default auxiliary convergence threshold.
    if isfinite(estimateRmse)
        threshold = max(1, 2 * estimateRmse);
    else
        threshold = 1;
    end

end

function timeVector = resolveTimeVector(numSamples, filterOutput, options)

    % Use the user-supplied time vector when it is available.
    if ~isempty(options.timeVector)
        timeVector = options.timeVector(:);
    % Otherwise build the time vector from the user-supplied step time.
    elseif ~isempty(options.stepTime)
        timeVector = (0:numSamples - 1).' * options.stepTime;
    % Otherwise build the time vector from the step time stored in the filter output.
    elseif isfield(filterOutput, 'stepTime')
        timeVector = (0:numSamples - 1).' * filterOutput.stepTime;
    else
        % Fall back to a sample index vector when no step time is available.
        timeVector = (0:numSamples - 1).';
    end

    % Check that the resolved time vector matches the number of samples.
    if numel(timeVector) ~= numSamples
        error('The resolved timeVector must match the number of samples.');
    end

end

function threshold = resolveConvergenceThreshold(fdTrue, fdMeas, options)

    % Use the user-supplied threshold when it is provided explicitly.
    if ~isempty(options.convergenceThresholdHz)
        threshold = options.convergenceThresholdHz;
        return
    end

    % Compute the measurement RMSE used to derive the default threshold.
    measurementRmse = trackingRmse(fdTrue, fdMeas);
    % Convert the RMSE into a default convergence threshold.
    if isfinite(measurementRmse)
        threshold = max(1, 2 * measurementRmse);
    else
        threshold = 1;
    end

end

