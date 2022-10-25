close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile  
[y,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_JONNE4_REAL.wav'); 

data = y(100000:end,1);

% Parameters
c = 299792458;                % Speed of light [m/s]
f_center = 5.8e9;             % Center Frequency [Hz]

% Processing
range = data*c/(4*pi*f_center); 

dataFFT = fft(normalize(data),size(data,1)*4);
dataFFT = dataFFT(1:end/4);
dataFFT = 20*log10(abs(dataFFT));
delta_f = linspace(0, Fs/2, size(data, 1));

% Plotting
figure(3), clf();
subplot(1,2,1); plot(range); title("Data Captured (Phase shift)"); grid on;

subplot(1,2,2); plot(delta_f, dataFFT); grid on; title("FFT"); xlim([0, 2]); ylim([70 100]);