
function [doppler, timeArr] = dopplerSearch(filename, Fs, timeSignal, Tcoh, t)

fc = 2.304e8;        % Carrier frequency [Hz] (passed to replica)
[r_local, ~] = localReplica(Fs, fc);

numTimeBlocks = timeSignal/Tcoh;

doppler = zeros(numTimeBlocks,1);
timeArr = zeros(numTimeBlocks,1);
fd = 0;

for n = 1:numTimeBlocks
    disp(n);
    disp(numTimeBlocks);    

    tStart = (n-1)*Tcoh;
    tEnd     = tStart + Tcoh; % End time [ms]
    [r, ~, ~, ~] = readSignal(filename, Fs, tStart, tEnd);
    
    r = r(:);
    fd_limit = 10000;

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fd*t);
    
    while fd_limit > 100
        
        % Correlate with local replica
        rcorr = xcorr(r_wiped, r_local);
    
        % Detection metric: peak magnitude
        likelihood = max(abs(rcorr));
        
        fd_upper = fd + fd_limit;
        r_wiped_upper= r .* exp(-1j*2*pi*fd_upper*t);
        rcorr_upper= xcorr(r_wiped_upper, r_local);
        likelihood_upper= max(abs(rcorr_upper));

        fd_lower = fd - fd_limit;
        r_wiped_lower= r .* exp(-1j*2*pi*fd_lower*t);
        rcorr_lower= xcorr(r_wiped_lower, r_local);
        likelihood_lower= max(abs(rcorr_lower));

        if likelihood_lower < likelihood && likelihood_upper < likelihood
            fd_limit = fd_limit/2;
            continue
        elseif likelihood_lower > likelihood
            fd = fd_lower;
            r_wiped = r .* exp(-1j*2*pi*fd*t);
            continue
        elseif likelihood_upper > likelihood
            fd = fd_upper;
            r_wiped = r .* exp(-1j*2*pi*fd*t);
        end
    end

    doppler(n) = fd;
    timeArr(n) = tStart;

end


return