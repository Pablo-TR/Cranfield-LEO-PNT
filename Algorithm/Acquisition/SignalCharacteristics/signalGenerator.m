%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            signalGenerator.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Algorithm
% Author: Alfredo
% Date: 8th March 2026
% Last Update: 8th March 2026
% Version: 1.0.0
%
% Description: Creation of synthetic signal, tReplicaStarts [ms]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [syntheticSignal, tReplicaStarts, time] = signalGenerator(Fs, timeSignal, fDoppler, SNR)
    
    timeSignal = timeSignal*1e-3;
    Fsym = 2.304e8;
    tRep = 1e-3; % Repeats every 1ms
    
    % r_local = localReplicaForGen(Fs, Fsym, SNR);
    [~, r_local] = localReplica(Fs, Fsym);

    r_local = r_local(:);     % ensure column
    L = length(r_local);              % replica length
    

    % timeLocal = (0:1/Fs:L/Fs-1/Fs)';
    % new_local = r_local .*exp(1j*2*pi.*fDoppler.*timeLocal);



    

    %

    samplesPerMs = round(Fs * tRep);          % 1 ms signal length
    
    % Create zero signal
    sig1ms = complex(zeros(samplesPerMs,1));


    % sig1ms = awgn(sig1ms, SNR);


    
    startSample = randi([1, samplesPerMs - L + 1]); % Random start position where replica fits
    
    tReplicaStarts = (startSample/Fs - 1/Fs)*1e3; % Time between 0 and 1ms where replica starts, -1/Fs because index 1 is for time 0secs and not 1/Fs
    
    % Insert replica
    sig1ms(startSample:startSample+L-1) = r_local;
    
    % Repeat signal for however much time specified (timeSignal)
    syntheticSignal = repmat(sig1ms, timeSignal*1e3, 1);

    % Constant Doppler addition
    time = (0:1/Fs:timeSignal-1/Fs)';
    syntheticSignal = syntheticSignal .*exp(1j*2*pi.*fDoppler.*time);

    % Noise addition
    syntheticSignal = awgn(syntheticSignal, SNR);
    
    %figure
    %plot(time, real(syntheticSignal))