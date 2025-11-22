% --- Initial Setup ---
clear all;
close all; % Close all previous figures
clc;       % Clear the command window

% --- File and Signal Parameters ---
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
Fs = 3.2e8; % Sampling frequency: 320 MHz
t = 10;    % Duration to read in milliseconds (ms)

% Calculate the total number of I/Q samples to read
L = floor(t * Fs / 1000);
N = L; % N is the number of I/Q samples

% --- File Read Parameters ---
data_type = 'int16';
bytes_per_sample = 2;

% --- Binary File Read ---
fprintf('Opening file: %s\n', filename);
[fid, msg] = fopen(filename, 'r');
if fid < 0
    error('Could not open file "%s" - Error: %s', filename, msg);
end

fprintf('Reading %d samples (%d ms)...\n', N, t);
data = fread(fid, [2, N], data_type);
fclose(fid);
fprintf('File read complete.\n');

if size(data, 2) < N
    warning('Fewer samples read than expected. End of file?');
    N = size(data, 2);
    L = N;
    if N == 0, error('No data read.'); end
end

% --- Signal Processing (Time Domain) ---
data_I = double(data(1, :)) / (2^15 - 1);
data_Q = double(data(2, :)) / (2^15 - 1);
signal = data_I + 1i * data_Q;
signal_detrended = detrend(signal);
fprintf('Signal processed (detrend).\n');

% --- FFT Calculation ---
signal_fft = fft(signal_detrended);
f_shifted = (-L/2 : L/2-1) * (Fs / L) * 1e-6;
fprintf('FFT calculation complete.\n');

% --- !!! CORRECTION FOR FIGURE 2 !!! ---
% --- PSD Calculation (pwelch) with AVERAGING ---
fprintf('Calculating PSD (pwelch) with averaging...\n');
% Define a SMALL window (segment) size.
% This will average many segments, reducing noise.
win_len = 8192; % 8192 samples per segment (you can adjust this)
window = hamming(win_len);
noverlap = win_len / 2; % 50% overlap
nfft = win_len; % Use the same number of FFT points as the window
% 'centered' gives us the spectrum from -Fs/2 to +Fs/2
% This function will now average (N / win_len) * 2 segments.
[pxx, f_pwelch] = pwelch(signal_detrended, window, noverlap, nfft, Fs, 'centered');
w = abs(f_pwelch(1)-f_pwelch(2));
f_pwelch = f_pwelch * 1e-6; % Convert frequency to MHz
fprintf('PSD (pwelch) calculation complete.\n');

% --- Main Plots ---
% % Plot 1: Centered FFT (Will look "noisy", this is normal)
% figure(1);
% plot(f_shifted, 10*log10(abs(fftshift(signal_fft))));
% title(sprintf('"Noisy" FFT (Raw) (%.1f ms of signal)', t));
% xlabel('Frequency [MHz]');
% ylabel('Amplitude [dB]');
% grid on;
% axis tight;

% Plot 2: PSD with pwelch (Will look "smooth", thanks to averaging)
figure(2);
plot(f_pwelch, 10*log10(pxx.*w)); 
title('"Smooth" PSD (pwelch with averaging)');
xlabel('Frequency [MHz]');
ylabel('Power/Frequency [dB/Hz]');
grid on;
axis tight;

% % --- Plot 3: 3D Spectrogram (smoothed version) ---
% fprintf('Calculating data for 3D spectrogram...\n');
% 
% % Parameters (longer window and overlap for smoother time/frequency transitions)
% win_spec = hamming(8192);    % larger window → smoother freq resolution
% overlap_spec = 0.75 * length(win_spec);  % higher overlap → smoother in time
% nfft_spec = 8192;
% 
% % Compute spectrogram
% [s_spec, f_spec, t_spec] = spectrogram(signal, win_spec, overlap_spec, nfft_spec, Fs, 'centered');
% t_spec_ms = t_spec * 1000;
% f_spec_mhz = f_spec * 1e-6;
% 
% % Convert to dB scale
% s_spec_db = 10 * log10(abs(s_spec));
% 
% % Optional numerical smoothing (moving average in time-frequency domain)
% s_spec_db_smooth = smoothdata(s_spec_db, 2, 'movmean', 5);  % smooth across time
% s_spec_db_smooth = smoothdata(s_spec_db_smooth, 1, 'movmean', 5); % smooth across frequency
% 
% % Plot 3D smooth spectrogram
% figure(3);
% surf(t_spec_ms, f_spec_mhz, s_spec_db_smooth, 'EdgeColor', 'none');
% shading interp;
% colormap jet;
% box on;
% axis tight;
% view(120, 30);
% title('3D Spectrogram (Smoothed)');
% xlabel('Time [ms]');
% ylabel('Frequency [MHz]');
% zlabel('Power [dB]');
% colorbar;
% 
% fprintf('Displaying 3D spectrogram.\n');
% 
% % --- Plot 4: Cyclic Autocorrelation (The "fingerprint") ---
% % --- THIS IS THE IMPROVED VERSION ---
% fprintf('Calculating cyclic autocorrelation (on a short chunk)...\n');
% % --- !!! HERE IS THE CHANGE !!! ---
% % Instead of using 'signal' (100 ms), we use a shorter chunk.
% % We will use only the first 10 ms (one OneWeb frame).
% t_chunk = 0.5; % 10 ms
% L_chunk = floor(t_chunk * Fs / 1000);
% signal_chunk = signal(1:L_chunk); % We take only the first chunk
% % The frequency range (alpha) remains the same
% alpha_hat = (230.39: 0.000001: 230.41) .* 1e6;
% % We calculate the autocorrelation ONLY on the short chunk
% [R, alpha] = cyclicAutoCorr(signal_chunk, 0, Fs, alpha_hat);
% % --- !!! NORMALIZATION !!! ---
% % Convert to magnitude and normalize (divide by the maximum)
% R_normalized = abs(R) / max(abs(R));
% % --- Plotting ---
% figure(4);
% plot(alpha_hat .* 1e-6, R_normalized, 'LineWidth', 1.5); % X-axis in MHz
% title('NORMALIZED Cyclic Spectrum (over 10ms)');
% xlabel('Cyclic Frequency (alpha) [MHz]');
% ylabel('Normalized Amplitude');
% grid on;
% axis tight;
% ylim([0, 1.1]); % Force Y-axis from 0 to ~1
% fprintf('Plot 4 (improved) complete.\n');
% 
% 
% figure(5);
% 
% plot(data_I,t);
% hold on
% plot(data_Q,t);
% hold off

