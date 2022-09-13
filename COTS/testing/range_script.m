close all; clear; clc;

%% COMPUTE RANGE FROM AUDIO FILE
% Read the audiofile 
[y,Fs] = audioread('Range_Test_File.m4a'); 

% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = y(:,2); % Sync data (square waveform)

% Parameters
c = 299792458;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
dr = c / (2*bandwidth);       % Range resolution [m]
Tp = 20e-3;                   % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
max_range = (N * dr)/2;       % Maximum range [m]

% Parse up-chirp data according to the sync data
sync_pulse = zeros(length(sync),2);
sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1
sync_pulse(:,2) = 1:1:length(sync_pulse(:,1)); % Set indexes

time_temp = find(sync_pulse==1);

% Build the matrix
up_data_parsed = zeros([], N); % Pre-allocate 
time = zeros(1,[]);            % Same
k = 1;
for i = 2:(size(sync_pulse)-N) 
    if sync_pulse(i,1) == 1 && sync_pulse(i-1) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = data(i:i+N-1)';
        time(1,k) = sync_pulse(i,2) / Fs;
        k = k + 1;
    end
end

% MS Clutter Rejection
final_data = bsxfun(@minus, up_data_parsed, mean(up_data_parsed, 1)); % Subtract column mean to each column

% 2-step MTI
MTI2 = zeros(size(up_data_parsed));
for t = 2:size(up_data_parsed,1)
    MTI2(t,:) = up_data_parsed(t,:) - up_data_parsed(t-1,:);
end

% 3-step MTI
MTI3 = zeros(size(up_data_parsed));
for t = 3:size(up_data_parsed,1)
    MTI3(t,:) = up_data_parsed(t,:) - 2*up_data_parsed(t-1,:) + up_data_parsed(t-2,:);
end

% FFT - No MS
sfft = 20*log10(abs(fft(up_data_parsed, 4*N, 2))); % dft using zero padding
sfft = sfft(:, 1:end/2) - max(max(sfft));          % Normalize data
      
% FFT - MS
sfft_ms = 20*log10(abs(fft(final_data, 4*N, 2)));  % dft using zero padding
sfft_ms = sfft_ms(:, 1:end/2) - max(max(sfft_ms)); % Normalize data      

% Plotting
f_start_plot = f_start/1e9;
f_stop_plot = f_stop/1e9;
range = linspace(0, max_range, 4*N);

%% Differnent BWs
f_start1 = 2.391e9;              % Start Frequency [Hz]
f_stop1 = 2.495e9;               % Stop Frequency  [Hz]
bandwidth1 = f_stop1 - f_start1; % [Hz]
dr1 = c / (2*bandwidth1);        % Range resolution [m]
max_range1 = (N * dr1)/2;        % Maximum range [m]
f_start_plot1 = f_start1/1e9;
f_stop_plot1 = f_stop1/1e9;
range1 = linspace(0, max_range1, 4*N);

f_start2 = 2.425e9;              % Start Frequency [Hz]
f_stop2 = 2.495e9;               % Stop Frequency  [Hz]
bandwidth2 = f_stop2 - f_start2; % [Hz]
dr2 = c / (2*bandwidth2);        % Range resolution [m]
max_range2 = (N * dr2)/2;        % Maximum range [m]
f_start_plot2 = f_start2/1e9;
f_stop_plot2 = f_stop2/1e9;
range2 = linspace(0, max_range2, 4*N);

figure;

a1=subplot(1,3,1);
imagesc(a1, range1, time, sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Larger BW, Tp="+Tp+"ms, f_{start}="+f_start_plot1+"GHz, f_{stop}="+f_stop_plot1+"GHz");

a2=subplot(1,3,2);
imagesc(a2, range, time, sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Normal BW, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a3=subplot(1,3,3);
imagesc(a3, range2, time, sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Narrower BW, Tp="+Tp+"ms, f_{start}="+f_start_plot2+"GHz, f_{stop}="+f_stop_plot2+"GHz");

%% MS - NO MS
figure;
a1=subplot(1,2,1);
imagesc(a1, range, time, sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-30 0]);
title("RTI without clutter rejection, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a2=subplot(1,2,2);
imagesc(a2, range, time, sfft_ms);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("RTI with clutter rejection, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

%% Different paddings
sfft1 = 20*log10(abs(fft(final_data, 2*N, 2))); 
sfft1 = sfft1(:, 1:end/2) - max(max(sfft1));          
range1 = linspace(0, max_range, 2*N);

sfft2 = 20*log10(abs(fft(final_data, 6*N, 2))); 
sfft2 = sfft2(:, 1:end/2) - max(max(sfft2));          
range2 = linspace(0, max_range, 6*N);

figure;
a1=subplot(1,3,1);
imagesc(a1, range1, time, sfft1);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Padding=2N, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a2=subplot(1,3,2);
imagesc(a2, range, time, sfft_ms);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Padding=4N, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a3=subplot(1,3,3);
imagesc(a3, range2, time, sfft2);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Padding=6N, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

%% MTI
% FFT - MTI2
sfft_MTI2 = 20*log10(abs(fft(MTI2, 4*N, 2)));  % dft using zero padding
sfft_MTI2 = sfft_MTI2(:, 1:end/2) - max(max(sfft_MTI2)); % Normalize data  

% FFT - MTI3
sfft_MTI3 = 20*log10(abs(fft(MTI3, 4*N, 2)));  % dft using zero padding
sfft_MTI3 = sfft_MTI3(:, 1:end/2) - max(max(sfft_MTI3)); % Normalize data  

figure;
a0=subplot(1,3,1);
imagesc(a0, range, time, sfft_ms);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("MS, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a1=subplot(1,3,2);
imagesc(a1, range, time, sfft_MTI2);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-40 0]);
title("2 point MTI, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a2=subplot(1,3,3);
imagesc(a2, range, time, sfft_MTI3);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-40 0]);
title("3 point MTI, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");