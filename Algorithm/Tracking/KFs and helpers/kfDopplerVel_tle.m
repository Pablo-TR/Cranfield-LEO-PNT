%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       kfDopplerVel.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Runs a 2-state Kalman filter with state [fd; fd_dot] and
% measurement vector [fd; fd_dot], where the rate measurement is derived
% from finite differences of the Doppler measurements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = kfDopplerVel_tle(fdMeas, fdDotModel, stepTime, sigmaFd, sigmaFddot,...
    processNoiseScale, measurementVar)

    % Read the total number of Doppler epochs.
    numEpochs = numel(fdMeas);
    % Define the 2-state constant-velocity transition matrix.
    Fmat = [1, stepTime;
            0, 1];
    % Define the 2-component measurement model for Doppler and modelled Doppler rate.
    Hmat = eye(2);

    % Build the finite-difference Doppler-rate sequence for diagnostics.
    fdDotMeas = zeros(numEpochs, 1);
    fdDotMeas(2:end) = diff(fdMeas) / stepTime;

    % Build the full measurement covariance for the coupled Doppler/rate measurement.
    R = [measurementVar^2, measurementVar^2 / stepTime;
         measurementVar^2 / stepTime, 2 * measurementVar^2 / stepTime^2];

    % Build the initial 2-state estimate.
    x = [fdMeas(1);
         fdDotMeas(2)];
    % Build the initial covariance matrix from the user-provided standard deviations.
    P = diag([(sigmaFd)^2, sigmaFddot^2]);

    % Build the process-noise matrix for the constant-velocity model.
    Q = processNoiseScale * [stepTime^3 / 3, stepTime^2 / 2;
                             stepTime^2 / 2, stepTime];

    % Preallocate the state history array.
    stateHistory = zeros(2, numEpochs);
    % Preallocate the raw innovation history for both measurement channels.
    innovation = zeros(2, numEpochs - 1);
    % Preallocate the normalised innovation history for both measurement channels.
    normalisedInnovation = zeros(2, numEpochs - 1);
    % Preallocate the covariance history array.
    covarianceHistory = zeros(2, 2, numEpochs);

    % Store the initial state estimate.
    stateHistory(:, 1) = x;
    % Store the initial covariance estimate.
    covarianceHistory(:, :, 1) = P;

    % Run the Kalman filter over the full Doppler sequence.
    for epoch = 1:numEpochs - 1
        xPred = Fmat * x;
        PPred = Fmat * P * Fmat' + Q;

        y = [fdMeas(epoch + 1);
             fdDotModel(epoch + 1)];
        S = Hmat * PPred * Hmat' + R;
        Kk = PPred * Hmat' / S;
        z = y - Hmat * xPred;
        cholFactor = chol(S, 'lower');
        nu = cholFactor \ z;

        x = xPred + Kk * z;
        josephTerm = eye(2) - Kk * Hmat;
        P = josephTerm * PPred * josephTerm' + Kk * R * Kk';

        stateHistory(:, epoch + 1) = x;
        innovation(:, epoch) = z;
        normalisedInnovation(:, epoch) = nu;
        covarianceHistory(:, :, epoch + 1) = P;
    end

    % Package the filter estimates and diagnostics into one output structure.
    output = struct( ...
        'fdEst', stateHistory(1, :).', ...
        'fdDotEst', stateHistory(2, :).', ...
        'fdDotMeas', fdDotMeas, ...
        'stateHistory', stateHistory, ...
        'innovation', innovation, ...
        'normalisedInnovation', normalisedInnovation, ...
        'covarianceHistory', covarianceHistory, ...
        'measurementCovariance', R, ...
        'stepTime', stepTime, ...
        'sigmaFd', sigmaFd, ...
        'processNoiseScale', processNoiseScale);

end

