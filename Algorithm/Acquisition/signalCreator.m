%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            signalCreator.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Algorithm
% Author: Alfredo
% Date: 8th March 2026
% Last Update: 8th March 2026
% Version: 1.0.0
%
% Description: Creation of synthetic signal
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [syntheticSignal, tReplicaStarts] = signalCreator(Fs, timeSignal, fDoppler, SNR)

    Fsym = 2.304e8;
    tRep = 1e-3; % Repeats every 1ms
    
    [r_local, ~] = localReplica(Fs, Fsym);
    
    r_local = r_local(:);     % ensure column
    L = length(r_local);              % replica length
    
    samplesPerMs = round(Fs * tRep);          % 1 ms signal length
    
    % Create zero signal
    sig1ms = complex(zeros(samplesPerMs,1));
    
    startSample = randi([1, samplesPerMs - L + 1]); % Random start position where replica fits
    
    tReplicaStarts = startSample/Fs - 1/Fs; % Time between 0 and 1ms where replica starts, -1/Fs because index 1 is for time 0secs and not 1/Fs
    
    % Insert replica
    sig1ms(startSample:startSample+L-1) = r_local;
    
    % Repeat signal for however much time specified (timeSignal)
    syntheticSignal = repmat(sig1ms, timeSignal, 1);

    % Constant Doppler addition
    time = (0:1/Fs:timeSignal*1e-3-1/Fs)';
    syntheticSignal = syntheticSignal .*exp(1j*2*pi.*fDoppler.*time);

    % Noise addition
    %noise = noiseStd * (randn(size(syntheticSignal)) + 1i*randn(size(syntheticSignal)));
    %syntheticSignal = syntheticSignal + noise;
    syntheticSignal = awgn(syntheticSignal, SNR);
    
    figure
    plot(time, real(syntheticSignal))