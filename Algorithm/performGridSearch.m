function [tau_found, fd_found] = performGridSearch(x_l, x_r, Fs)
fd_max = 614*10^3;
tau_max = 4e-3;
fd = -fd_max:3070:fd_max;
tau = 
% Circshift (account for time delay) of local replica

% Perform wipe-off to received signal
    x_wiped = x_r .* exp(-1j*2*pi*fd.*k.*Ts);
% Compute FFT of wiped received signal
    x_wiped_fft = fft(x_wiped, K);
% Cross-correlate with time-shifted local replica
    Y = x_wiped_fft * x_local_fft_shifted;
% Compute ambiguity function
    R = (1/K^2) * ifft(Y, K);
% Obtain maximum and indices
    S  = max(R.^2);
    [i, j] = find(S == R.^2);
    fd_found = fd(i);
    tau_found = tau(j);
% Plot ambiguity function
figure
surf(fd, tau, R.^2)
end