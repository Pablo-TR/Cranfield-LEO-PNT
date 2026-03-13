function localExtended = localExtendedReplica(Fs, timeSignal)


    Fsym = 2.304e8;
    tRep = 1e-3; % Repeats every 1ms
    
    [~, r_local] = localReplica(Fs, Fsym);
    
    r_local = r_local(:);     % ensure column
    L = length(r_local);              % replica length
    
    samplesPerMs = round(Fs * tRep);          % 1 ms signal length
    
    % Create zero signal
    sig1ms = complex(zeros(samplesPerMs,1));

    % Insert replica
    sig1ms(1:L) = r_local;
    
    % Repeat signal for however much time specified (timeSignal)
    localExtended = repmat(sig1ms, timeSignal, 1);

    time = (0:1/Fs:timeSignal*1e-3-1/Fs)';

    

    figure
    plot(time, real(localExtended))