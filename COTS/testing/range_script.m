close all; clear all;

%% COMPUTE RANGE FROM AUDIO FILE
% Read the audiofile 
[y,Fs] = audioread('Range_Test_File.m4a'); 

% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

% Parameters
c = 299792458;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
Tp = 20e-3;                   % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
dr = c / (2*bandwidth);       % Range resolution [m]
max_range = (N * dr)/2;       % Maximum range [m]

% Parse up-chirp data according to the sync data
sync_pulse = zeros(length(sync),2);
sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1
sync_pulse(:,2) = 1:1:length(sync_pulse(:,1)); % Set indexes

time_temp = find(sync_pulse==1);

% Vectorized approach -> not working
% up_data = data(sync_pulse(:,1)==1);         % Use sync_pulse as a mask on the original data to get the up-chirp values
% 
% up_data_first = up_data(1:N:end);
% % up_data_first(:,2) = time_temp(1:N:end);
% ind_list = time_temp(1:N:end);
% pulse_data = [];
% for i=2:length(up_data_first)
%     if ind_list(i)+N-1<=N
%         first_el = up_data_first(i);
%         chirp = [first_el, data(ind_list(i):ind_list(i)+N-1)'];
%         pulse_data = cat(1, pulse_data, chirp);
%     end
% end
% 
% X = mod(-mod(length(pulse_data), N), N);         % Used to find the previous divisible value with respect to length(up_data)
% up_data_cut = pulse_data(1:end-(N-X));           % Remove the last elements so that we can reshape up_data
% up_data_parsed = reshape(up_data_cut, [], N)';
% 
% time = zeros(1,length(up_data_parsed));
% j = 1;
% for i = 1:N+1:length(time_temp)
%     time(j) = time_temp(i) / Fs;
%     j = j+1;
% end

% Classic for loop
up_data_parsed = zeros([], N); % Pre-allocate matrix
time = zeros(1,[]);            % Same
k = 1;
for i = 11:(size(sync_pulse)-N) 
    if sync_pulse(i,1) == 1 && mean(sync_pulse(i-10:i-1,1)) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = data(i:i+N-1)';
        time(1,k) = sync_pulse(i,2) / Fs;
        k = k + 1;
    end
end

% MS Clutter Rejection
final_data = bsxfun(@minus, up_data_parsed, mean(up_data_parsed, 1)); % Subtract column mean to each column

% FFT - No MS
sfft = 20*log10(abs(fft(up_data_parsed, 4*N, 2))); % dft using zero padding
sfft = sfft(:, 1:end/2) - max(max(sfft));          % Normalize data
      
% FFT - MS
sfft_ms = 20*log10(abs(fft(final_data, 4*N, 2)));  % dft using zero padding
sfft_ms = sfft_ms(:, 1:end/2) - max(max(sfft_ms)); % Normalize data  

range = linspace(0, max_range, 4*N);    

% Plotting
f_start_plot = f_start/1e9;
f_stop_plot = f_stop/1e9;

a1=subplot(1,2,1);
imagesc(a1, range, time, sfft);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-50 0]);
title("RTI without clutter rejection, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

a2=subplot(1,2,2);
imagesc(a2, range, time, sfft_ms);
xlim([0 100]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-50 0]);
title("RTI with clutter rejection, Tp="+Tp+"ms, f_{start}="+f_start_plot+"GHz, f_{stop}="+f_stop_plot+"GHz");

