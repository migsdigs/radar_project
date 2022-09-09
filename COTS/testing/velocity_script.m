close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
[y,Fs] = audioread('Velocity_Test_File.m4a'); 

% Take care of data inversion by the sound card
data = -y(:,1); % Intensity of the received signal

% Parameters
c = 299792458;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

% Parse the data
X = mod(-mod(length(data), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data(1:end-(N-X));           % Remove the last elements so that we can reshape up_data
data_parsed = reshape(data_cut, [], N);

% MS Clutter Rejection
final_data = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract row mean to each column

% FFT across the rows for multiple pulses
% 4*N padding for smoother visualization
f = abs(fft(final_data, 4*N, 2));
f = 20*log10(f);

% Only considering lower half of the frequency domain
f = f(:,1:size(f, 2) / 2);
f = f - max(max(f));
delta_f = linspace(0, Fs/2, size(f,2));
vel = (delta_f*c)/(2 * f_center);
time = linspace(1, Tp * size(f,1), size(f,1));

% Visualization of Velocity information
imagesc(vel,time,f);
caxis([-35 0]);
colorbar;
set(gca,'XLim',[0 40]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
toc