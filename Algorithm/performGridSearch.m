function [tau_found, fd_found] = performGridSearch(rxl, rxr, Fs)
N = size(rx) 
k = round(tau*Fs); % Compute samples for circular shift (sample delay)
 rxl_shift = circshift(rxl, k); % Shift local replica by tau seconds
 rxl_shift_doppler = rxl_shift * exp(1j*fd.*tvec); 
end