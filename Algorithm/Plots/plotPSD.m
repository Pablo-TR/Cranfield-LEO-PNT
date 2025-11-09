function plotPSD(pxx, fpsd)
    plot(fpsd, 10*log10(pxx)); 
    title('"Smooth" PSD (pwelch with averaging)');
    xlabel('Frequency [MHz]');
    ylabel('Power/Frequency [dB/Hz]');
    grid on;
    axis tight;
end