function plotPSD(fpsd, pxx)
    plot(fpsd, 10*log10(pxx)); 
    title('"Smooth" PSD (pwelch with averaging)');
    xlabel('Frequency [MHz]');
    ylabel('Power [dB]');
    grid on;
    axis tight;
end