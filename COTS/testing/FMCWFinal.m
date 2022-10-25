clear, clc
% Read the audiofile
[y,Fs] = audioread('audacity_recordings\single_target_range.wav');%single_target_range.wav'); 
numRidges = 1; %Number of targets

% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

% Parameters
c = 299792458;                % Speed of light [m/s]
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
f_center = (f_start + f_stop)/2;
Tp = 20e-3;                   % Pulse width (up-chirp) [s]
N = Tp * Fs;                  % Number of samples per up-chrip pulse
dr = c / (2*bandwidth);       % Range resolution [m]
max_range = (N * dr)/2;       % Maximum range [m]

fridgeLimitRange = 25;
fridgeLengthRange = round((2*bandwidth*fridgeLimitRange*4)/c); % To limit the search of tfridge
freqsRange = linspace(fridgeLimitRange, 0, fridgeLengthRange);

fridgeLimitVel = 10;
fridgeLengthVel = round((2*f_center*fridgeLimitVel*4)/c); % To limit the search of tfridge
freqsVel = linspace(fridgeLimitVel, 0, fridgeLengthVel);

range = linspace(0, max_range, 4*N);

% Parse up-chirp data according to the sync data
sync_pulse = zeros(length(sync),2);
sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1
sync_pulse(:,2) = 1:1:length(sync_pulse(:,1)); % Set indexes

time_temp = find(sync_pulse==1);

% Build the matrix where timesteps are row-wise
up_data_parsed = zeros([], N); % Pre-allocate 
down_data_parsed = zeros([],N);
time = zeros(1,[]);            % Same
k = 1;
for i = 2:(size(sync_pulse)-2*N) 
    if sync_pulse(i,1) == 1 && sync_pulse(i-1) == 0 % First value of a row = first up-chirp value
        up_data_parsed(k,:) = data(i:i+N-1)';
        down_data_parsed(k,:) = data(i+N:i+2*N-1)';
        time(1,k) = sync_pulse(i,2) / Fs;
        k = k + 1;
    end
end

% MS Clutter Rejection
up_data_parsed = bsxfun(@minus, up_data_parsed, mean(up_data_parsed, 1)); % Subtract column mean to each column
down_data_parsed = bsxfun(@minus, down_data_parsed, mean(down_data_parsed, 1)); % Subtract column mean to each column

% 3-step MTI
MTI3_up = zeros(size(up_data_parsed));
for t = 3:size(up_data_parsed,1)
    MTI3_up(t,:) = up_data_parsed(t,:) - 2*up_data_parsed(t-1,:) + up_data_parsed(t-2,:);
end

MTI3_down = zeros(size(down_data_parsed));
for t = 3:size(down_data_parsed,1)
    MTI3_down(t,:) = down_data_parsed(t,:) - 2*down_data_parsed(t-1,:) + down_data_parsed(t-2,:);
end

% FFT
sfft_up = 20*log10(abs(fft(MTI3_up, 4*N, 2)));  % dft using zero padding
sfft_up = sfft_up(:, 1:end/2) - max(max(sfft_up)); % Normalize data 
sfft_up(sfft_up < -1000000) = -1000000; % Make sure no -Inf values messing it up

% FFT
sfft_down = 20*log10(abs(fft(MTI3_down, 4*N, 2)));  % dft using zero padding
sfft_down = sfft_down(:, 1:end/2) - max(max(sfft_down)); % Normalize data 
sfft_down(sfft_down < -1000000) = -1000000; % Make sure no -Inf values messing it up

sfftFridgeUp = fft(MTI3_up, 4*N,2);
sfftFridgeDown = fft(MTI3_down, 4*N,2);

% find fridges
[fridge_up, ~, ~] = tfridge(rot90(sfftFridgeUp(:, 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',numRidges);
% [fridge1, ~, ~] = tfridge(rot90(sfftFridgeUp(1:floor(end*6.5/13), 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',numRidges,'NumFrequencyBins',15);
% [fridge2, ~, ~] = tfridge(rot90(sfftFridgeUp((floor(end*6.5/13)+1):floor(end*9.3/13), 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',numRidges,'NumFrequencyBins',1);
% [fridge3, ~, ~] = tfridge(rot90(sfftFridgeUp((floor(end*9.3/13)+1):end, 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',1);
% fridge_up = [fridge1; fridge2; [fridge3, zeros(length(fridge3),1)]];
[fridge_down, ~, ~] = tfridge(rot90(sfftFridgeDown(:, 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',numRidges); %,'NumFrequencyBins',1);
% [fridge11, ~, ~] = tfridge(rot90(sfftFridgeDown(1:floor(end*6.5/13), 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',numRidges,'NumFrequencyBins',15);
% [fridge22, ~, ~] = tfridge(rot90(sfftFridgeDown((floor(end*6.5/13)+1):floor(end*9.3/13), 1:fridgeLengthRange)), freqsRange, 2,'NumRidges',numRidges,'NumFrequencyBins',1);
% [fridge33, ~, ~] = tfridge(rot90(sfftFridgeDown((floor(end*9.3/13)+1):end, 1:fridgeLengthRange)), freqsRange, 1,'NumRidges',1);
% fridge_down = [fridge11; fridge22; [fridge33, zeros(length(fridge33),1)]];

% compute range and velocity
f_beat = (fridge_down + fridge_up) / 2;
rangeExtracted = 15*(f_beat * c * Tp) / bandwidth;
f_dopler = (fridge_down - fridge_up) / 2;
velExtracted = 20*(c * f_dopler) / (2 * f_center);

% Compute range and velocity using difference method
% rangeExtracted = 15*(fridge_up*c*Tp)/bandwidth;
% velExtracted = zeros(numRidges,length(rangeExtracted));
% velExtracted(:,5:end) = (rangeExtracted(5:end,:)' - rangeExtracted(1:end-4,:)')/(5*Tp);

%% Plotting
% Assembling the figure
figure();
sgtitle('Single target simultaneous range and velocity measurement using difference method') 
spectogram=subplot(1,3,1);
imagesc(range, time, sfft_up);
xlim([0 25]);
xlabel('Range (m)');
ylabel('Time (s)');
colorbar;
caxis([-35 0]);
title("Range-spectogram");

subplot(1,3,2)
rangePlot = plot(time,rangeExtracted);
title("Range vs Time plot"); xlabel('Time [s]'); ylabel('Range [m]');

subplot(1,3,3)
velocityPlot = plot(time,velExtracted);
title("Velocity vs Time plot"); xlabel('Time [s]'); ylabel('velocity [m/s]');
