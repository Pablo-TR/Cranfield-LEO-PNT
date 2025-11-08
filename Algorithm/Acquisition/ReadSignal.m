% --- Configuración Inicial ---
clear all;
close all; % Cierra todas las figuras anteriores
clc;       % Limpia la ventana de comandos

% --- Parámetros de Archivo y Señal ---
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
Fs = 3.2e8; % Frecuencia de muestreo: 320 MHz
t = 1;      % Duración a leer en milisegundos (ms)

% Calcular el número total de muestras I/Q a leer
L = floor(t * Fs / 1000);
N = L; % N es el número de muestras I/Q

% --- Parámetros de Lectura de Archivo ---
% MEJORA: Asunción crítica sobre el tipo de datos.
% La mayoría de los SDR guardan como 'int16' (entero 16-bit) o 'single'.
% El 'uint8' por defecto de fread() es casi seguro incorrecto.
% Probaremos con 'int16'.
data_type = 'int16';
bytes_per_sample = 2; % 'int16' usa 2 bytes

% --- Lectura del Archivo Binario ---
[fid, msg] = fopen(filename, 'r');
if fid < 0
    error('No se pudo abrir el archivo "%s" - Error: %s', filename, msg);
end

% CORRECCIÓN: Leer los datos especificando el tipo 'int16'.
% El formato [2, N] lee [I, Q, I, Q, ...] y los organiza en 2 filas
data = fread(fid, [2, N], data_type);
fclose(fid);

% MEJORA: Comprobar si se leyeron suficientes datos
if size(data, 2) < N
    warning('Se leyeron menos muestras de las esperadas (N=%d). Fin de archivo?', size(data, 2));
    N = size(data, 2); % Actualizar N al número real de muestras leídas
    if N == 0
        error('No se leyeron datos.');
    end
    L = N;
end

% --- Procesamiento de Señal (Dominio del Tiempo) ---
% Convertir a double y normalizar (asumiendo 16-bit con signo)
data_I = double(data(1, :)) / (2^15 - 1);
data_Q = double(data(2, :)) / (2^15 - 1);
signal = data_I + 1i * data_Q;

% CORRECCIÓN: Eliminar el offset de DC en el dominio del tiempo.
% Esto es lo que causa el gran pico en 0 Hz ("DC Spike").
% Usar detrend() es la forma correcta, no eliminar el bin de la FFT.
signal_detrended = detrend(signal);

% --- Cálculo de FFT ---
% Usar la señal sin offset ('signal_detrended')
signal_fft = fft(signal_detrended);

% CORRECCIÓN: Crear los vectores de frecuencia correctos.
% Vector de frecuencia para FFT centrada (-Fs/2 a +Fs/2) - en MHz
f_shifted = (-L/2 : L/2-1) * (Fs / L) * 1e-6;

% --- Cálculo de PSD (pwelch) ---
% MEJORA: pwelch es más robusto para estimar la PSD.
% Lo centramos para que coincida con la FFT.
window = hamming(N);
noverlap = N / 2; % Solapamiento del 50%
nfft = max(256, 2^nextpow2(N)); % Tamaño de FFT

% 'centered' nos da directamente el espectro de -Fs/2 a +Fs/2
[pxx, f_pwelch] = pwelch(signal_detrended, window, noverlap, nfft, Fs, 'centered');

% Convertir frecuencia a MHz
f_pwelch = f_pwelch * 1e-6;

% --- Gráficas ---

% Gráfica 1: FFT Centrada (La más útil)
figure(1);
% MEJORA: Usar fftshift() y el vector f_shifted.
% MEJORA: Plotear en dB (10*log10) es estándar para espectros.
plot(f_shifted, 10*log10(abs(fftshift(signal_fft))));
title(sprintf('FFT Centrada (%.1f ms de señal)', t));
xlabel('Frecuencia [MHz]');
ylabel('Amplitud [dB]');
grid on;
axis tight; % Ajustar ejes

% Gráfica 2: PSD con pwelch (Más suave y precisa)
figure(2);
% CORRECCIÓN: No multiplicar pxx por w. Solo convertir a dB.
plot(f_pwelch, 10*log10(pxx)); 
title('Power Spectral Density (pwelch)');
xlabel('Frecuencia [MHz]');
ylabel('Potencia/Frecuencia [dB/Hz]');
grid on;
axis tight;

% Gráfica 3: Espectrograma (Muy útil para ver la señal en el tiempo)
% MEJORA: Una FFT de 1ms es un solo "snapshot". Un espectrograma
% te muestra si la señal está cambiando o saltando.
figure(3);
% Usamos ventanas más pequeñas para ver la evolución en el tiempo
win_spec = hamming(4096);
overlap_spec = 2048;
nfft_spec = 4096;
spectrogram(signal, win_spec, overlap_spec, nfft_spec, Fs, 'centered', 'yaxis');
title('Espectrograma');
xlabel('Tiempo');
% El 'yaxis' ya pone la frecuencia en MHz o GHz automáticamente
ylabel('Frecuencia');


% --- Autocorrelación Cíclica (Opcional) ---
% Esta sección requiere una función 'cyclicAutoCorr' que no es estándar de MATLAB.
% Si tienes esta función en tu "path", descomenta las siguientes líneas.

% disp('Calculando autocorrelación cíclica...');
% alpha_hat = (230.39: 0.0005: 230.41) .* 1e6;
% [R, alpha] = cyclicAutoCorr(signal, 0, Fs, alpha_hat);
% 
% figure(4); % Crear una nueva figura
% plot(alpha_hat .* 1e-6, abs(R));
% title('Espectro Cíclico (Autocorrelación)');
% xlabel('Frecuencia Cíclica (alpha) [MHz]');
% ylabel('Amplitud');
% grid on;
