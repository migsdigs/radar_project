clear, clc

% Read the audiofile
[y,Fs] = audioread('audacity_recordings\single_target_range.wav'); 
numRidges = 1; % Amout of targets

% In case of IQ Demodulation, choose which data to use from here
% [I,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_REAL.wav'); 
% [Q,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_IMAG.wav'); 
% % Composite the complex value and its conjugate
% data = complex(I,Q);
% data = conj(data);

% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

% Parameters
c = 299792458;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
dr = c / (2*bandwidth);       % Range resolution [m]
Tp = 20e-3;                   % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
max_range = (N * dr)/2;       % Maximum range [m]

fridgeLimit = 25;
fridgeLength = round((2*bandwidth*fridgeLimit*4)/c); % To limit the search of tfridge
freqs = linspace(fridgeLimit, 0, fridgeLength);
%freqs = linspace(Fs / (4 / 2), 0, 4*N); % To limit the search of tfridge
range = linspace(0, max_range, 4*N);

% Parse up-chirp data according to the sync data
sync_pulse = zeros(length(sync),2);
sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1
sync_pulse(:,2) = 1:1:length(sync_pulse(:,1)); % Set indexes

time_temp = find(sync_pulse==1);

% Build the matrix where timesteps are row-wise
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

% 3-step MTI
MTI3 = zeros(size(up_data_parsed));
for t = 3:size(up_data_parsed,1)
    MTI3(t,:) = up_data_parsed(t,:) - 2*up_data_parsed(t-1,:) + up_data_parsed(t-2,:);
end

% FFT - MTI3
sfft_MTI3 = 20*log10(abs(fft(MTI3, 4*N, 2)));  % dft using zero padding
sfft_MTI3 = sfft_MTI3(:, 1:end/2) - max(max(sfft_MTI3)); % Normalize data 
sfft_MTI3(sfft_MTI3 < -1000000) = -1000000; % Make sure no -Inf values messing it up

% find fridges
sfft_fridge = fft(MTI3, 4*N, 2);
[fridge, ~, ~] = tfridge(rot90(sfft_fridge(:, 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',10);
% [fridge1, ~, ~] = tfridge(rot90(sfft_fridge(1:floor(end*7.2/13), 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',10);
% [fridge2, ~, ~] = tfridge(rot90(sfft_fridge((floor(end*7.2/13)+1):floor(end*9.3/13), 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',1);
% [fridge3, ~, ~] = tfridge(rot90(sfft_fridge((floor(end*9.3/13)+1):end, 1:fridgeLength)), freqs, 1,'NumRidges',1);
% fridge = [fridge1; fridge2; [fridge3, zeros(length(fridge3),1)]];
% Convert to range
rangeExtracted = 15*(fridge * c * Tp) / bandwidth;

% Assembling the figure
figure();
sgtitle('Single target range measurement') 
spectogram=subplot(1,2,1);
imagesc(range, time, sfft_MTI3);
xlim([0 25]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Range-spectogram");

subplot(1,2,2)
rangePlot = plot(time,rangeExtracted);
title("Range vs Time plot"); xlabel('Time [s]'); ylabel('Range [m]');
