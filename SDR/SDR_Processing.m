close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
% [y,Fs] = audioread('Velocity_Test_File.m4a'); 
[I,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_REAL.wav'); 
[Q,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_IMAG.wav'); 
nTargets = 1;

% Composite the complex value and its conjugate
data1 = complex(I,Q);
data2 = conj(data1);

% Parameters
c = 299792458;                % Speed of light [m/s]
f_center = 5.8e9;             % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

% Parse the data
X = mod(-mod(length(data1), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data1((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';
final_data1 = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract column mean to each column

X = mod(-mod(length(data2), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data2((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';
final_data2 = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % MS Clutter rejection

% FFT 
f1 = abs(fft(final_data1, 4*N, 2));
f1 = 20*log10(f1);
f1 = f1(:,1:size(f1, 2) / 2);

f2 = abs(fft(final_data2, 4*N, 2));
f2 = 20*log10(f2);
f2 = f2(:,1:size(f2, 2) / 2);

delta_f = linspace(0, Fs/2, size(f1, 2)); 
vel = (delta_f * c)/(2 * f_center);
time = linspace(1, Tp * size(f1, 1), size(f1, 1));

% Normalization
% f1_norm = f1 - max(max(f1));
% f2_norm = f2 - max(max(f2));

f_row_max = max(f1, [], 2);
f1_norm = f1 - f_row_max;

f_row_max = max(f2, [], 2);
f2_norm = f2 - f_row_max;

% find fridges
[fridge1, ~, ~] = tfridge(rot90(f1(:, 1:82)), delta_f(1:82), 1,'NumRidges',nTargets);
[fridge2, ~, ~] = tfridge(rot90(f2(:, 1:82)), delta_f(1:82), 1,'NumRidges',nTargets);
vel1 = (c * fridge1) / (2 * f_center);
vel2 = (c * fridge2) / (2 * f_center);

% Plot
figure(1);
a1=subplot(2,2,1);
imagesc(vel, time, f1_norm);
caxis([-10 0]);
colorbar;
set(gca,'XLim',[0 7]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("I + jQ");

a2=subplot(2,2,2);
imagesc(vel, time, f2_norm);
caxis([-10 0]);
colorbar;
set(gca,'XLim',[0 7]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("I - jQ");

a3 = subplot(2,2,3);
plot(time, vel1);

a4 = subplot(2,2,4);
plot(time, vel2);

