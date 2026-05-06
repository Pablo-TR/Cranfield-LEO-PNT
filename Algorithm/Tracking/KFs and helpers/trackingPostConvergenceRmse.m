%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               trackingPostConvergenceRmse.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 15th March 2026
% Last Update: 15th March 2026
% Version: 1.0.0
%
% Description: Computes the RMSE from the first convergence instant to the
% end of the run, using the same convergence definition as
% trackingConvergenceTime.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metric = trackingPostConvergenceRmse(reference, estimate, ...
    timeVector, threshold, sustainSamples)

    % Use five consecutive samples as the default convergence requirement.
    if nargin < 5
        sustainSamples = 5;
    end

    % Reuse the convergence metric so the same convergence criterion is applied.
    convergenceMetric = trackingConvergenceTime(reference, estimate, ...
        timeVector, threshold, sustainSamples);

    % Initialise the output structure with the convergence summary.
    metric = struct( ...
        'threshold', threshold, ...
        'sustainSamples', sustainSamples, ...
        'hasConverged', convergenceMetric.hasConverged, ...
        'convergenceIndex', convergenceMetric.convergenceIndex, ...
        'convergenceTime', convergenceMetric.convergenceTime, ...
        'sampleCount', 0, ...
        'rmseValue', NaN);

    % Exit early when the sequence never converges.
    if ~convergenceMetric.hasConverged
        return
    end

    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);

    % Extract the reference tail from the convergence instant onward.
    postConvReference = reference(convergenceMetric.convergenceIndex:end);
    % Extract the estimate tail from the convergence instant onward.
    postConvEstimate = estimate(convergenceMetric.convergenceIndex:end);
    % Keep only the finite sample pairs inside the post-convergence window.
    validMask = isfinite(postConvReference) & isfinite(postConvEstimate);

    % Store the number of valid post-convergence samples.
    metric.sampleCount = nnz(validMask);
    % Exit early when no valid post-convergence sample remains.
    if metric.sampleCount == 0
        return
    end

    % Compute the RMSE over the post-convergence valid samples.
    metric.rmseValue = trackingRmse(postConvReference(validMask), ...
        postConvEstimate(validMask));

end

