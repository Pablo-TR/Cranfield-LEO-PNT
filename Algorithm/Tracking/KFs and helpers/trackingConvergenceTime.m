%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 trackingConvergenceTime.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Computes the first time at which the tracking error remains
% inside a specified tolerance for a sustained number of samples.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metric = trackingConvergenceTime(reference, estimate, ...
    timeVector, threshold, sustainSamples)

    % Use five consecutive samples as the default convergence requirement.
    if nargin < 5
        sustainSamples = 5;
    end

    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);
    % Convert the time axis into a column vector.
    timeVector = timeVector(:);

    % Initialise the output structure with the default non-converged state.
    metric = struct( ...
        'threshold', threshold, ...
        'sustainSamples', sustainSamples, ...
        'hasConverged', false, ...
        'convergenceIndex', NaN, ...
        'convergenceTime', NaN);

    % Compute the absolute tracking error at every sample.
    errorAbs = abs(estimate - reference);
    % Mark the samples where the absolute error is finite.
    validMask = isfinite(errorAbs);
    % Mark the samples that satisfy the convergence threshold.
    isWithinThreshold = validMask & (errorAbs <= threshold);

    % Exit early when the sequence is shorter than the sustain window.
    if numel(reference) < sustainSamples
        return
    end

    % Search for the first sample from which the error remains inside the threshold.
    for sampleIdx = 1:numel(reference) - sustainSamples + 1
        if all(isWithinThreshold(sampleIdx:sampleIdx + sustainSamples - 1))
            % Mark the sequence as converged.
            metric.hasConverged = true;
            % Store the first sustained-convergence sample index.
            metric.convergenceIndex = sampleIdx;
            % Store the corresponding convergence time.
            metric.convergenceTime = timeVector(sampleIdx);
            return
        end
    end

end

