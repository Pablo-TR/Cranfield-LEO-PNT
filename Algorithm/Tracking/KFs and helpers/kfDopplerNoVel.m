%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      kfDopplerNoVel.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Runs a 2-state Kalman filter with state [fd; fd_dot] and a
% scalar Doppler measurement.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = kfDopplerNoVel(fdMeas, stepTime, FdSigma, FddotSigma,...
    processNoiseScale, measurementVar, max_acc)
    % Read the total number of Doppler epochs.
    numEpochs = numel(fdMeas);
    % Define the 2-state constant-velocity transition matrix.
    Fmat = [1, stepTime;
            0, 1];
    % Define the scalar Doppler measurement model.
    Hmat = [1, 0];

    % Initialise the Doppler-rate state from the first finite difference when possible.
    if numEpochs >= 2
        initialFdDot = (fdMeas(2) - fdMeas(1)) / stepTime;
    else
        % Fall back to zero rate when only one sample is available.
        initialFdDot = 0;
    end

    % Build the initial state vector.
    x = [fdMeas(1);
         initialFdDot];

    % Use the covariance model selected by the caller.
    if nargin < 7 || isempty(max_acc)
        % Build the initial covariance from the user-provided standard deviations.
        P = diag([FdSigma^2, FddotSigma^2]);
        % Build the continuous-white-acceleration process-noise matrix.
        Q = processNoiseScale * [stepTime^3 / 3, stepTime^2 / 2;
                                 stepTime^2 / 2, stepTime];
    else
        % Convert the maximum acceleration prior into state-scale bounds.
        s1 = 0.5 * max_acc * stepTime^2;
        s2 = max_acc * stepTime;
        % Build an alternative initial covariance using the measurement scale.
        P = diag([measurementVar^2, max(initialFdDot^2, eps)]);
        % Build the diagonal process-noise matrix from the acceleration bounds.
        Q = diag([s1^2, s2^2]);
    end

    % Preallocate the state history array.
    stateHistory = zeros(2, numEpochs);
    % Preallocate the raw innovation history.
    innovation = zeros(1, numEpochs - 1);
    % Preallocate the normalised innovation history.
    normalisedInnovation = zeros(1, numEpochs - 1);
    % Preallocate the covariance history array.
    covarianceHistory = zeros(2, 2, numEpochs);

    % Store the initial state estimate.
    stateHistory(:, 1) = x;
    % Store the initial covariance estimate.
    covarianceHistory(:, :, 1) = P;

    % Build the scalar measurement-noise variance.
    R = measurementVar^2;

    % Run the Kalman filter over the full Doppler sequence.
    for epoch = 1:numEpochs - 1
        xPred = Fmat * x;
        PPred = Fmat * P * Fmat' + Q;

        y = fdMeas(epoch + 1);
        S = Hmat * PPred * Hmat' + R;
        Kk = (PPred * Hmat') / S;
        z = y - Hmat * xPred;
        nu = z / sqrt(S);

        x = xPred + Kk * z;
        josephTerm = eye(2) - Kk * Hmat;
        P = josephTerm * PPred * josephTerm' + Kk * R * Kk';

        stateHistory(:, epoch + 1) = x;
        innovation(epoch) = z;
        normalisedInnovation(epoch) = nu;
        covarianceHistory(:, :, epoch + 1) = P;
    end

    % Package the filter estimates and diagnostics into one output structure.
    output = struct( ...
        'fdEst', stateHistory(1, :).', ...
        'fdDotEst', stateHistory(2, :).', ...
        'stateHistory', stateHistory, ...
        'innovation', innovation(:), ...
        'normalisedInnovation', normalisedInnovation(:), ...
        'covarianceHistory', covarianceHistory, ...
        'stepTime', stepTime, ...
        'measurementVar', measurementVar^2, ...
        'processNoiseScale', processNoiseScale);

end

