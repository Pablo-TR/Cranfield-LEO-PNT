%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               mainTracking2DTests.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 14th April 2026
% Last Update: 15th April 2026
% Version: 1.1.0
%
% Description: Main file to test the 2D Doppler model and the 2-state and
% 3-state linear Kalman filters in a controlled environment.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);
set(groot, 'defaultLegendFontSize', 14);
set(groot, 'defaultLineLineWidth', 1.5);

% Fix the random generator for repeatable results.
rng(7);

% Define the baseline model timestep in seconds.
baseStepTime = 1e-3;

% Define the orbital height in metres.
baseH = 1200e3;

% Define the carrier frequency in Hz.
baseF0 = 11.325e9;

% Define the mixing frequency in Hz.
baseFmix = -1;

% Define the initial central angle in rad.
baseLambda0 = -0.6;

% Define the baseline measurement noise standard deviation in Hz.
baseSigmaNoise = 1000;

% Define the baseline process-noise scale.
baseProcessNoiseScale = 10;

% Define the Doppler convergence threshold in Hz.
dopplerConvergenceThreshold = 10;

% Define the number of consecutive samples required for convergence.
consecutiveSamplesRequired = 5;

% Generate the baseline 2D Doppler model.
[tBase, fdBase, dfdtBase, d2fdt2Base] = generateOWDoppler2DModel(baseStepTime, 0, baseH, baseF0, baseFmix, baseLambda0);

% Convert the time vector to a column vector.
tBase = tBase(:);

% Build a Doppler-rate ground-truth vector with the same length as fdBase.
fdDotBase = [dfdtBase(1); dfdtBase(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST 1. Validation of filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the ground truth as the measurement vector.

fdMeasNominal = addNoise(fdBase, baseSigmaNoise);

% Run the 2-state Kalman filter.
output2Nominal = linearKF(fdMeasNominal, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

% Extract the estimated Doppler from the 2-state filter.
fdEst2Nominal = output2Nominal.fdEst;

% Extract the estimated Doppler rate from the 2-state filter.
fdDotEst2Nominal = output2Nominal.fdDotEst;

% Run the 3-state Kalman filter.
output3Nominal = linearKF(fdMeasNominal, baseStepTime, baseProcessNoiseScale,baseSigmaNoise, d2fdt2Base, 1);

% Extract the estimated Doppler from the 3-state filter.
fdEst3Nominal = output3Nominal.fdEst;

% Extract the estimated Doppler rate from the 3-state filter.
fdDotEst3Nominal = output3Nominal.fdDotEst;

% Compute the Doppler error of the 2-state filter.
fdError2Nominal = fdEst2Nominal - fdBase;

% Compute the Doppler-rate error of the 2-state filter.
fdDotError2Nominal = fdDotEst2Nominal - fdDotBase;

% Compute the Doppler error of the 3-state filter.
fdError3Nominal = fdEst3Nominal - fdBase;

% Compute the Doppler-rate error of the 3-state filter.
fdDotError3Nominal = fdDotEst3Nominal - fdDotBase;

% Compute the Doppler RMSE of the 2-state filter.
rmseFd2Nominal = sqrt(mean(fdError2Nominal.^2));

% Compute the Doppler-rate RMSE of the 2-state filter.
rmseFdDot2Nominal = sqrt(mean(fdDotError2Nominal.^2));

% Compute the Doppler RMSE of the 3-state filter.
rmseFd3Nominal = sqrt(mean(fdError3Nominal.^2));

% Compute the Doppler-rate RMSE of the 3-state filter.
rmseFdDot3Nominal = sqrt(mean(fdDotError3Nominal.^2));

% Initialise the convergence index of the 2-state filter.
convIndex2Nominal = NaN;

% Search for the first epoch from which the Doppler error remains within the convergence band.
for k = 1:(length(fdBase) - consecutiveSamplesRequired + 1)
    fdErrorWindow = abs(fdEst2Nominal(k:k + consecutiveSamplesRequired - 1) - fdBase(k:k + consecutiveSamplesRequired - 1));
    if all(fdErrorWindow <= dopplerConvergenceThreshold)
        convIndex2Nominal = k;
        break
    end
end

% Initialise the convergence index of the 3-state filter.
convIndex3Nominal = NaN;

% Search for the first epoch from which the Doppler error remains within the convergence band.
for k = 1:(length(fdBase) - consecutiveSamplesRequired + 1)
    fdErrorWindow = abs(fdEst3Nominal(k:k + consecutiveSamplesRequired - 1) - fdBase(k:k + consecutiveSamplesRequired - 1));
    if all(fdErrorWindow <= dopplerConvergenceThreshold)
        convIndex3Nominal = k;
        break
    end
end

% Convert the 2-state convergence index into time.
if isnan(convIndex2Nominal)
    convTime2Nominal = NaN;
else
    convTime2Nominal = tBase(convIndex2Nominal);
end

% Convert the 3-state convergence index into time.
if isnan(convIndex3Nominal)
    convTime3Nominal = NaN;
else
    convTime3Nominal = tBase(convIndex3Nominal);
end

% Display a header for the first test.
disp(' ')
disp('TEST 1. NOMINAL CASE')

% Display the 2-state metrics.
disp(['2-state KF Doppler RMSE [Hz]: ', num2str(rmseFd2Nominal)])
disp(['2-state KF Doppler-rate RMSE [Hz/s]: ', num2str(rmseFdDot2Nominal)])
disp(['2-state KF convergence time [s]: ', num2str(convTime2Nominal)])

% Display the 3-state metrics.
disp(['3-state KF Doppler RMSE [Hz]: ', num2str(rmseFd3Nominal)])
disp(['3-state KF Doppler-rate RMSE [Hz/s]: ', num2str(rmseFdDot3Nominal)])
disp(['3-state KF convergence time [s]: ', num2str(convTime3Nominal)])

% Create a figure for the nominal Doppler estimation.
figure;
plot(tBase, fdBase.*1e-3, 'k', 'LineWidth', 1.5)
hold on
plot(tBase, fdEst2Nominal.*1e-3, '--', 'LineWidth', 1.5)
plot(tBase, fdEst3Nominal.*1e-3, '-.', 'LineWidth', 2)
plot(tBase, fdMeasNominal.*1e-3, '--', 'LineWidth', 2)
grid on
xlabel('Time [s]')
ylabel('Doppler [Hz]')
title('2nd order KF validation')
legend('Ground truth', '2-state KF', '3-state KF','Measurements', 'Location', 'best')


figure;
plot(tBase, fdDotBase, 'k', 'LineWidth', 1.5)
hold on
plot(tBase, fdDotEst3Nominal, '-.', 'LineWidth', 1)
plot(tBase, fdDotEst2Nominal, '--', 'LineWidth', 2)
grid on
xlabel('Time [s]')
ylabel('Doppler rate [Hz/s]')
title('Doppler-rate estimation')
legend('Ground truth', '2-state KF', '3-state KF','Measurements', 'Location', 'best')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4. COMPARISON BETWEEN REGIONS OF THE S-CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the representative centre times in seconds.
regionCenterTimeVec = [400 650 850];

% Define the region labels.
regionLabels = {'Early', 'Middle', 'Late'};

% Define the fixed observation-window duration in seconds.
regionWindowDuration = 100;

% Set the number of Monte Carlo runs used in this test.
test4NumRuns = 100;

% Store the number of analysed regions.
numRegions = length(regionCenterTimeVec);

% Preallocate the region start indices.
regionStartIndexVec = zeros(numRegions, 1);

% Preallocate the region end indices.
regionEndIndexVec = zeros(numRegions, 1);

% Preallocate the local time vectors for each region.
tWindowCell = cell(numRegions, 1);

% Preallocate the local time vectors starting at zero for each region.
tWindowLocalCell = cell(numRegions, 1);

% Preallocate the Doppler truth windows for each region.
fdWindowCell = cell(numRegions, 1);

% Preallocate the Doppler-rate truth windows for each region.
fdDotWindowCell = cell(numRegions, 1);

% Precompute the indices and truth windows of each region.
for i = 1:numRegions

    % Find the index nearest to the desired centre time.
    [~, centreIndex] = min(abs(tBase - regionCenterTimeVec(i)));

    % Compute half of the window length in samples.
    halfWindowSamples = round(regionWindowDuration / (2 * baseStepTime));

    % Compute the first index of the region.
    startIndex = centreIndex - halfWindowSamples;

    % Compute the last index of the region.
    endIndex = centreIndex + halfWindowSamples;

    % Limit the first index to the available signal.
    startIndex = max(startIndex, 1);

    % Limit the last index to the available signal.
    endIndex = min(endIndex, length(fdBase));

    % Store the first index of the region.
    regionStartIndexVec(i) = startIndex;

    % Store the last index of the region.
    regionEndIndexVec(i) = endIndex;

    % Extract the local time vector.
    tWindowCell{i} = tBase(startIndex:endIndex);

    % Shift the local time to start at zero.
    tWindowLocalCell{i} = tWindowCell{i} - tWindowCell{i}(1);

    % Extract the ground-truth Doppler in the region.
    fdWindowCell{i} = fdBase(startIndex:endIndex);

    % Extract the ground-truth Doppler rate in the region.
    fdDotWindowCell{i} = fdDotBase(startIndex:endIndex);
end

% Preallocate the 2-state Doppler RMSE by region and run.
rmseFd2RegionRuns = zeros(numRegions, test4NumRuns);

% Preallocate the 2-state Doppler-rate RMSE by region and run.
rmseFdDot2RegionRuns = zeros(numRegions, test4NumRuns);

% Preallocate the 2-state convergence time by region and run.
convTime2RegionRuns = NaN(numRegions, test4NumRuns);

% Preallocate the 3-state Doppler RMSE by region and run.
rmseFd3RegionRuns = zeros(numRegions, test4NumRuns);

% Preallocate the 3-state Doppler-rate RMSE by region and run.
rmseFdDot3RegionRuns = zeros(numRegions, test4NumRuns);

% Preallocate the 3-state convergence time by region and run.
convTime3RegionRuns = NaN(numRegions, test4NumRuns);

% Preallocate the example 2-state Doppler estimates for each region.
exampleFdEst2Cell = cell(numRegions, 1);

% Preallocate the example 3-state Doppler estimates for each region.
exampleFdEst3Cell = cell(numRegions, 1);

% Run the Monte Carlo analysis for the regions of the S-curve.
for runIdx = 1:test4NumRuns

    % Build one noisy measurement vector for the current Monte Carlo run.
    fdMeasRegion = addNoise(fdBase, baseSigmaNoise);

    % Loop over the different regions.
    for i = 1:numRegions

        % Extract the first index of the current region.
        startIndex = regionStartIndexVec(i);

        % Extract the last index of the current region.
        endIndex = regionEndIndexVec(i);

        % Extract the local time vector.
        tWindow = tWindowCell{i};

        % Extract the local time vector starting at zero.
        tWindowLocal = tWindowLocalCell{i};

        % Extract the ground-truth Doppler in the region.
        fdWindow = fdWindowCell{i};

        % Extract the noisy Doppler measurements in the region.
        fdMeasWindow = fdMeasRegion(startIndex:endIndex);

        % Extract the ground-truth Doppler rate in the region.
        fdDotWindow = fdDotWindowCell{i};

        % Run the 2-state Kalman filter on the region.
        output2 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

        % Extract the 2-state Doppler estimate.
        fdEst2 = output2.fdEst;

        % Extract the 2-state Doppler-rate estimate.
        fdDotEst2 = output2.fdDotEst;

        % Run the 3-state Kalman filter on the region.
        output3 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

        % Extract the 3-state Doppler estimate.
        fdEst3 = output3.fdEst;

        % Extract the 3-state Doppler-rate estimate.
        fdDotEst3 = output3.fdDotEst;

        % Compute the 2-state Doppler RMSE for the current region and run.
        rmseFd2RegionRuns(i, runIdx) = sqrt(mean((fdEst2 - fdWindow).^2));

        % Compute the 2-state Doppler-rate RMSE for the current region and run.
        rmseFdDot2RegionRuns(i, runIdx) = sqrt(mean((fdDotEst2 - fdDotWindow).^2));

        % Compute the 3-state Doppler RMSE for the current region and run.
        rmseFd3RegionRuns(i, runIdx) = sqrt(mean((fdEst3 - fdWindow).^2));

        % Compute the 3-state Doppler-rate RMSE for the current region and run.
        rmseFdDot3RegionRuns(i, runIdx) = sqrt(mean((fdDotEst3 - fdDotWindow).^2));

        % Search for 2-state convergence within the region.
        for k = 1:(length(fdWindow) - consecutiveSamplesRequired + 1)
            fdErrorWindow = abs(fdEst2(k:k + consecutiveSamplesRequired - 1) - fdWindow(k:k + consecutiveSamplesRequired - 1));
            if all(fdErrorWindow <= dopplerConvergenceThreshold)
                convTime2RegionRuns(i, runIdx) = tWindowLocal(k);
                break
            end
        end

        % Search for 3-state convergence within the region.
        for k = 1:(length(fdWindow) - consecutiveSamplesRequired + 1)
            fdErrorWindow = abs(fdEst3(k:k + consecutiveSamplesRequired - 1) - fdWindow(k:k + consecutiveSamplesRequired - 1));
            if all(fdErrorWindow <= dopplerConvergenceThreshold)
                convTime3RegionRuns(i, runIdx) = tWindowLocal(k);
                break
            end
        end

        % Keep the first Monte Carlo run as the example run for plotting.
        if runIdx == 1
            exampleFdEst2Cell{i} = fdEst2;
            exampleFdEst3Cell{i} = fdEst3;
        end
    end
end

% Preallocate the mean 2-state Doppler RMSE by region.
rmseFd2Region = zeros(numRegions, 1);

% Preallocate the mean 2-state Doppler-rate RMSE by region.
rmseFdDot2Region = zeros(numRegions, 1);

% Preallocate the mean 2-state convergence time by region.
convTime2Region = NaN(numRegions, 1);

% Preallocate the mean 3-state Doppler RMSE by region.
rmseFd3Region = zeros(numRegions, 1);

% Preallocate the mean 3-state Doppler-rate RMSE by region.
rmseFdDot3Region = zeros(numRegions, 1);

% Preallocate the mean 3-state convergence time by region.
convTime3Region = NaN(numRegions, 1);

% Compute the mean Monte Carlo metrics of each region.
for i = 1:numRegions

    % Compute the mean 2-state Doppler RMSE.
    rmseFd2Region(i) = mean(rmseFd2RegionRuns(i, :));

    % Compute the mean 2-state Doppler-rate RMSE.
    rmseFdDot2Region(i) = mean(rmseFdDot2RegionRuns(i, :));

    % Compute the mean 3-state Doppler RMSE.
    rmseFd3Region(i) = mean(rmseFd3RegionRuns(i, :));

    % Compute the mean 3-state Doppler-rate RMSE.
    rmseFdDot3Region(i) = mean(rmseFdDot3RegionRuns(i, :));

    % Extract the finite 2-state convergence times.
    validConvTime2 = convTime2RegionRuns(i, isfinite(convTime2RegionRuns(i, :)));

    % Compute the mean 2-state convergence time when available.
    if ~isempty(validConvTime2)
        convTime2Region(i) = mean(validConvTime2);
    end

    % Extract the finite 3-state convergence times.
    validConvTime3 = convTime3RegionRuns(i, isfinite(convTime3RegionRuns(i, :)));

    % Compute the mean 3-state convergence time when available.
    if ~isempty(validConvTime3)
        convTime3Region(i) = mean(validConvTime3);
    end
end

% Create a figure to show the analysed regions on the S-curve.
figure;

% Plot the complete Doppler curve.
plot(tBase, fdBase.*1e-3, 'k', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Loop over the different regions.
for i = 1:numRegions
    % Plot the analysed region on top of the complete S-curve.
    plot(tWindowCell{i}, fdWindowCell{i}.*1e-3, 'LineWidth', 1.5)
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [kHz]')

% Add title.
title('Test 4. Regions selected on the S-curve')

% Add legend.
legend('Complete S-curve', 'Early region', 'Middle region', 'Late region', 'Location', 'best')

% Create a figure for the region comparison metrics.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE by region.
bar(categorical(regionLabels), [rmseFd2Region rmseFd3Region].*1e-3)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [kHz]')

% Add title.
title(sprintf('Test 4. Mean comparison between regions of the S-curve (%d runs)', test4NumRuns))

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'northwest')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE by region.
bar(categorical(regionLabels), [rmseFdDot2Region rmseFdDot3Region].*1e-3)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [kHz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time by region.
bar(categorical(regionLabels), [convTime2Region convTime3Region])

% Add grid.
grid on

% Label the vertical axis.
ylabel('Convergence time [s]')

% Create a figure for one representative run with the 2-state KF.
figure;

% Plot the complete Doppler curve.
plot(tBase, fdBase.*1e-3, 'k', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Loop over the different regions.
for i = 1:numRegions

    % Plot the 2-state KF estimate of the current region.
    if i == 1
        plot(tWindowCell{i}, exampleFdEst2Cell{i}.*1e-3, '--', ...
            'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, ...
            'DisplayName', '2-state KF estimate')
    else
        plot(tWindowCell{i}, exampleFdEst2Cell{i}.*1e-3, '--', ...
            'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, ...
            'HandleVisibility', 'off')
    end
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [kHz]')

% Add title.
title('Test 4. Representative run with the 2-state KF')

% Add legend.
legend('Ground truth', '2-state KF estimate', 'Location', 'best')

% Create a figure for one representative run with the 3-state KF.
figure;

% Plot the complete Doppler curve.
plot(tBase, fdBase.*1e-3, 'k', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Loop over the different regions.
for i = 1:numRegions

    % Plot the 3-state KF estimate of the current region.
    if i == 1
        plot(tWindowCell{i}, exampleFdEst3Cell{i}.*1e-3, '-.', ...
            'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, ...
            'DisplayName', '3-state KF estimate')
    else
        plot(tWindowCell{i}, exampleFdEst3Cell{i}.*1e-3, '-.', ...
            'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, ...
            'HandleVisibility', 'off')
    end
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [kHz]')

% Add title.
title('Test 4. Representative run with the 3-state KF')

% Add legend.
legend('Ground truth', '3-state KF estimate', 'Location', 'best')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 6. SENSITIVITY TO THE SAMPLING PERIOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the list of sampling periods in seconds.
stepTimeVec = [1e-3 1e-2 1e-1 1 10];

% Preallocate the 2-state Doppler RMSE vector.
rmseFd2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 2-state Doppler-rate RMSE vector.
rmseFdDot2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 2-state convergence-time vector.
convTime2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state Doppler RMSE vector.
rmseFd3Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state Doppler-rate RMSE vector.
rmseFdDot3Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state convergence-time vector.
convTime3Step = zeros(length(stepTimeVec), 1);

% Loop over the different sampling periods.
for i = 1:length(stepTimeVec)

    % Extract the current sampling period.
    currentStepTime = stepTimeVec(i);

    % Generate the 2D Doppler model for the current sampling period.
    [tCurrent, fdCurrent, dfdtCurrent, d2fdt2Current] = generateOWDoppler2DModel(currentStepTime, 0, baseH, baseF0, baseFmix, baseLambda0);
    % Convert the current time vector to a column vector.
    tCurrent = tCurrent(:);

    % Build a Doppler-rate ground-truth vector with the same length as fdCurrent.
    fdDotCurrent = [dfdtCurrent(1); dfdtCurrent(:)];

    % Add Gaussian noise to the current measurements.
    fdMeasCurrent = addNoise(fdCurrent, baseSigmaNoise);

    % Run the 2-state Kalman filter for the current sampling period.
    output2 = linearKF(fdMeasCurrent, currentStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Current, 0);

    % Extract the 2-state Doppler estimate.
    fdEst2 = output2.fdEst;

    % Extract the 2-state Doppler-rate estimate.
    fdDotEst2 = output2.fdDotEst;

    % Run the 3-state Kalman filter for the current sampling period.
    output3 = linearKF(fdMeasCurrent, currentStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Current, 1);

    % Extract the 3-state Doppler estimate.
    fdEst3 = output3.fdEst;

    % Extract the 3-state Doppler-rate estimate.
    fdDotEst3 = output3.fdDotEst;

    % Compute the 2-state Doppler RMSE.
    rmseFd2Step(i) = sqrt(mean((fdEst2 - fdCurrent).^2));

    % Compute the 2-state Doppler-rate RMSE.
    rmseFdDot2Step(i) = sqrt(mean((fdDotEst2 - fdDotCurrent).^2));

    % Compute the 3-state Doppler RMSE.
    rmseFd3Step(i) = sqrt(mean((fdEst3 - fdCurrent).^2));

    % Compute the 3-state Doppler-rate RMSE.
    rmseFdDot3Step(i) = sqrt(mean((fdDotEst3 - fdDotCurrent).^2));

    % Initialise the current 2-state convergence time.
    convTime2Step(i) = NaN;

    % Search the 2-state convergence time.
    for k = 1:(length(fdCurrent) - consecutiveSamplesRequired + 1)
        fdErrorWindow = abs(fdEst2(k:k + consecutiveSamplesRequired - 1) - fdCurrent(k:k + consecutiveSamplesRequired - 1));
        if all(fdErrorWindow <= dopplerConvergenceThreshold)
            convTime2Step(i) = tCurrent(k);
            break
        end
    end

    % Initialise the current 3-state convergence time.
    convTime3Step(i) = NaN;

    % Search the 3-state convergence time.
    for k = 1:(length(fdCurrent) - consecutiveSamplesRequired + 1)
        fdErrorWindow = abs(fdEst3(k:k + consecutiveSamplesRequired - 1) - fdCurrent(k:k + consecutiveSamplesRequired - 1));
        if all(fdErrorWindow <= dopplerConvergenceThreshold)
            convTime3Step(i) = tCurrent(k);
            break
        end
    end
end

% Create a figure for the sampling-period sensitivity.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE versus sampling period.
semilogx(stepTimeVec, rmseFd2Step, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state Doppler RMSE versus sampling period.
semilogx(stepTimeVec, rmseFd3Step, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 6. Sensitivity to the sampling period')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'best')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE versus sampling period.
semilogx(stepTimeVec, rmseFdDot2Step, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state Doppler-rate RMSE versus sampling period.
semilogx(stepTimeVec, rmseFdDot3Step, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time versus sampling period.
semilogx(stepTimeVec, convTime2Step, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state convergence time versus sampling period.
semilogx(stepTimeVec, convTime3Step, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Sampling period [s]')

% Label the vertical axis.
ylabel('Convergence time [s]')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 8. SENSITIVITY TO THE FILTER TUNING USING ONLY THE Q/R FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the list of Q/R fractions to analyse.
qOverRVec = [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1 3 10 30 100 1000];

% Build one noisy measurement vector for this test.
fdMeasTuning = addNoise(fdBase, baseSigmaNoise);

% Preallocate the 2-state Doppler RMSE vector.
rmseFd2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 2-state Doppler-rate RMSE vector.
rmseFdDot2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 2-state convergence-time vector.
convTime2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state Doppler RMSE vector.
rmseFd3Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state Doppler-rate RMSE vector.
rmseFdDot3Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state convergence-time vector.
convTime3Tuning = zeros(length(qOverRVec), 1);

% Loop over the different Q/R fractions.
for i = 1:length(qOverRVec)

    % Extract the current Q/R ratio.
    currentQoverR = qOverRVec(i);

    % Convert the Q/R ratio into the process-noise scale expected by linearKF.
    currentProcessNoiseScale = currentQoverR * baseSigmaNoise^2;

    % Run the 2-state Kalman filter with the current Q/R ratio.
    output2 = linearKF(fdMeasTuning, baseStepTime, currentProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

    % Extract the 2-state Doppler estimate.
    fdEst2 = output2.fdEst;

    % Extract the 2-state Doppler-rate estimate.
    fdDotEst2 = output2.fdDotEst;

    % Run the 3-state Kalman filter with the current Q/R ratio.
    output3 = linearKF(fdMeasTuning, baseStepTime, currentProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

    % Extract the 3-state Doppler estimate.
    fdEst3 = output3.fdEst;

    % Extract the 3-state Doppler-rate estimate.
    fdDotEst3 = output3.fdDotEst;

    % Compute the 2-state Doppler RMSE.
    rmseFd2Tuning(i) = sqrt(mean((fdEst2 - fdBase).^2));

    % Compute the 2-state Doppler-rate RMSE.
    rmseFdDot2Tuning(i) = sqrt(mean((fdDotEst2 - fdDotBase).^2));

    % Compute the 3-state Doppler RMSE.
    rmseFd3Tuning(i) = sqrt(mean((fdEst3 - fdBase).^2));

    % Compute the 3-state Doppler-rate RMSE.
    rmseFdDot3Tuning(i) = sqrt(mean((fdDotEst3 - fdDotBase).^2));

    % Initialise the 2-state convergence time.
    convTime2Tuning(i) = NaN;

    % Search the 2-state convergence time.
    for k = 1:(length(fdBase) - consecutiveSamplesRequired + 1)
        fdErrorWindow = abs(fdEst2(k:k + consecutiveSamplesRequired - 1) - fdBase(k:k + consecutiveSamplesRequired - 1));
        if all(fdErrorWindow <= dopplerConvergenceThreshold)
            convTime2Tuning(i) = tBase(k);
            break
        end
    end

    % Initialise the 3-state convergence time.
    convTime3Tuning(i) = NaN;

    % Search the 3-state convergence time.
    for k = 1:(length(fdBase) - consecutiveSamplesRequired + 1)
        fdErrorWindow = abs(fdEst3(k:k + consecutiveSamplesRequired - 1) - fdBase(k:k + consecutiveSamplesRequired - 1));
        if all(fdErrorWindow <= dopplerConvergenceThreshold)
            convTime3Tuning(i) = tBase(k);
            break
        end
    end
end

% Create a figure for the tuning sensitivity.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE versus Q/R.
semilogx(qOverRVec, rmseFd2Tuning, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state Doppler RMSE versus Q/R.
semilogx(qOverRVec, rmseFd3Tuning, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 8. Sensitivity to the Q/R fraction')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'best')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE versus Q/R.
semilogx(qOverRVec, rmseFdDot2Tuning, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state Doppler-rate RMSE versus Q/R.
semilogx(qOverRVec, rmseFdDot3Tuning, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time versus Q/R.
semilogx(qOverRVec, convTime2Tuning, '--o', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 3-state convergence time versus Q/R.
semilogx(qOverRVec, convTime3Tuning, '-s', 'LineWidth', 1.5)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Q/R fraction')

% Label the vertical axis.
ylabel('Convergence time [s]')
