function plotSpectogram(t,f,s)
    surf(t, f, s);
    shading interp;
    box on;
    axis tight;
    view(120, 30);
    title('3D Spectrogram (Surf Plot)');
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    zlabel('Power [dB]');
    fprintf('Displaying plots.\n');
end
