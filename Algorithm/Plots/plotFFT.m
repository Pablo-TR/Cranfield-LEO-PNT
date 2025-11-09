function plotFFT(fShifted, signalFFT)
    figure(1);
    plot(fShifted, 10*log10(abs(fftshift(signalFFT))));
    title('"Noisy" FFT');
    xlabel('Frequency [MHz]');
    ylabel('Amplitude [dB]');
    grid on;
    axis tight;
end