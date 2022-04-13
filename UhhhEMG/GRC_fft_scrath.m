[pxx,f] = pspectrum(data.emg(16,:));
figure;
plot(f,pow2db(pxx))
grid on
xlabel('Frequency (Hz)')
ylabel('Power Spectrum (dB)')
title('overground')


%% via fft
Fs = 2000;
x =  data.emg(16,:);
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;
figure;
plot(freq,10*log10(psdx))
grid on
title('No Stim Walking: Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')