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
data_cut = data((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';

% MS Clutter Rejection
final_data = bsxfun(@minus, data_parsed, mean(data_parsed, 1)); % Subtract column mean to each column

% FFT 
f = abs(fft(final_data, 4*N, 2));
f = 20*log10(f);

f = f(:,1:size(f, 2) / 2);

% Normalization 1
f1 = f - max(max(f));
delta_f1 = linspace(0, Fs/2, size(f1,2)); 
vel1 = (delta_f1*c)/(2 * f_center);
time1 = linspace(1, Tp * size(f1,1), size(f1,1));

% Normalization 2
f2 = bsxfun(@minus, f, max(f, 2)); 
delta_f2 = linspace(0, Fs/2, size(f2,2)); 
vel2 = (delta_f2*c)/(2 * f_center);
time2 = linspace(1, Tp * size(f2,1), size(f2,1));

% Plotting
imagesc(vel1, time1, f1);
caxis([-35 0]);
colorbar;
set(gca,'XLim',[0 40]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
