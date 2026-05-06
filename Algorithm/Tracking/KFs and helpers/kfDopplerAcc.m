%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       kfDopplerAcc.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Runs a 3-state Kalman filter with state
% [fd; fd_dot; fd_ddot] and a scalar Doppler measurement.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = kfDopplerAcc(fdMeas, stepTime, FdSigma, FddotSigma,...
    processNoiseScale, measurementVar)

    % Read the total number of Doppler epochs.
    numEpochs = numel(fdMeas);
    % Define the 3-state constant-acceleration transition matrix.
    Fmat = [1, stepTime, 0.5 * stepTime^2;
            0, 1, stepTime;
            0, 0, 1];
    % Define the scalar Doppler measurement model.
    Hmat = [1, 0, 0];

    % Build the initial 3-state estimate from the first samples.
    x = [fdMeas(1);
         (fdMeas(2)-fdMeas(1))/stepTime;
         0];
    % Build the initial covariance matrix.
    P = diag([(FdSigma)^2, (FddotSigma)^2, (1e-2)^2]);

    % Build the process-noise matrix for the constant-acceleration model.
    Q = processNoiseScale * [stepTime^5 / 20, stepTime^4 / 8, ...
                             stepTime^3 / 6;
                             stepTime^4 / 8, stepTime^3 / 3, ...
                             stepTime^2 / 2;
                             stepTime^3 / 6, stepTime^2 / 2, ...
                             stepTime];

    % Preallocate the state history array.
    stateHistory = zeros(3, numEpochs);
    % Preallocate the raw innovation history.
    innovation = zeros(1, numEpochs - 1);
    % Preallocate the normalised innovation history.
    normalisedInnovation = zeros(1, numEpochs - 1);
    % Preallocate the covariance history array.
    covarianceHistory = zeros(3, 3, numEpochs);

    % Store the initial state estimate.
    stateHistory(:, 1) = x;
    % Store the initial covariance estimate.
    covarianceHistory(:, :, 1) = P;

    % Run the Kalman filter over the full Doppler sequence.
    for epoch = 1:numEpochs - 1
        xPred = Fmat * x;
        PPred = Fmat * P * Fmat' + Q;

        y = fdMeas(epoch + 1);
        S = Hmat * PPred * Hmat' + measurementVar^2;
        Kk = (PPred * Hmat') / S;
        z = y - Hmat * xPred;
        nu = z / sqrt(S);

        x = xPred + Kk * z;
        josephTerm = eye(3) - Kk * Hmat;
        P = josephTerm * PPred * josephTerm' + ...
            Kk * measurementVar^2 * Kk';

        stateHistory(:, epoch + 1) = x;
        innovation(epoch) = z;
        normalisedInnovation(epoch) = nu;
        covarianceHistory(:, :, epoch + 1) = P;
    end

    % Package the filter estimates and diagnostics into one output structure.
    output = struct( ...
        'fdEst', stateHistory(1, :).', ...
        'fdDotEst', stateHistory(2, :).', ...
        'fdDdotEst', stateHistory(3, :).', ...
        'stateHistory', stateHistory, ...
        'innovation', innovation(:), ...
        'normalisedInnovation', normalisedInnovation(:), ...
        'covarianceHistory', covarianceHistory, ...
        'stepTime', stepTime, ...
        'measurementVar', measurementVar, ...
        'processNoiseScale', processNoiseScale);

end

