clear all
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';   

%t0 = 10; % ms
%t1 = 100; %ms
t = 1; %ms

%bytes_per_sample = 8; % int16 = 2 bytes
Fs = 3.2e8;         % Sampling frequency: 320 Megasamples/s
%L1 = floor((t0/1000) * Fs); % Sample number where to start reading
%L2 = floor((t1/1000) * Fs);% Sample number where to end reading
%N = L2-L1 + 1; % Number of samples to be read
L = floor(t * Fs / 1000); % Total number of samples to read
N = L; % Number of samples to be read
[fid, msg] = fopen(filename, 'r');
%fseek(fid, (L1-1) * 2 * bytes_per_sample, 'bof');
%data = fread(fid, [2, N], 'double');
data = fread(fid, [2, N]);
fclose(fid);

if fid < 0
    error('Failed to open file "%s" because "%s"', filename, msg);
end



data_I = data(1, :); % I component
data_Q = data(2, :); % Q component 

%data_db = 10*log10(data_I);
%time = 0:(1/Fs):seconds_to_read;
%plot(time(1:end-1), real(data_db))

signal = data_I + data_Q.*1i;
signal_fft = fft(signal);
f = Fs/L*(0:L-1); 
f = f * 1e-6; % frequency vector in MHz.

%Removing data point f=0 MHz due to offset?
f = f(2:end);
signal_fft = signal_fft(2:end);


%% PSD computation
signal_detrended = detrend(signal);
window = hamming(length(signal_detrended));
signal_windowed = signal_detrended .* window';
[pxx,w] = pwelch(signal_windowed, window,[], [],Fs);
dw = w(2)-w(1);
pxx = pxx.*w;
%% Autocorrelation
alpha_hat = (230.39: 0.0005: 230.41).*1e6;
[R, alpha] = cyclicAutoCorr(signal,0,Fs,alpha_hat);

plot(alpha_hat,abs(R))

plotFFT(f, abs(signal_fft), time_to_read, pxx, w)


%% Plots
function plotFFT(f, signal_fft, t, pxx, w)
    figure(1)
    plot(f, abs(signal_fft));
    title(sprintf('FFT for %d miliseconds of signal', t))
    xlabel('|f| [MHz]')
    ylabel('Amplitude')
    
    figure(2)
    plot(f, abs(fftshift(signal_fft)));
    title(sprintf('FFT for %d miliseconds of signal', t))
    xlabel('|f| [MHz]')
    ylabel('Amplitude')

    figure(3)
    plot(w, abs(10*log10(pxx)));
    title('Power Spectral Density (PSD)');
    xlabel('Frequency [Hz]');
    ylabel('Power [dB]');
    grid on;
end

