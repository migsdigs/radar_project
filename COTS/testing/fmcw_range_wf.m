close all; clear; clc;

%% COMPUTE RANGE FROM AUDIO FILE
% Read the audiofile 
[y,Fs] = audioread('audacity_recordings\single_target_range.wav'); 

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
Tp = 20e-3;                   % Pulse width (up-chirp) [s]
N = Tp * Fs;                  % Number of samples per up-chrip pulse
max_range = (N * dr)/2;       % Maximum range [m]

wf = 128;

%% 

% % plotting
spec_range = ones(wf, 2*N) * -Inf;

range = linspace(0, max_range, N*2);
data_t = linspace(0, wf*Tp*2, wf);

% create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
spectrogram = imagesc(range, data_t, spec_range, [-20 0]); colorbar; set(gca,'XLim',[0 100]);
title("Range vs Time Spectrogram"); xlabel('Velocity [m/s]'); ylabel('Time [s]');


index = 1;
start_index = 1;
end_index = 0;
upchrip = false;
pulse = zeros(1,[]);
x = true;

while x == true
    if sync(index+1,1) > 0 && sync(index,1) < 0
        start_index = index + 1;
        end_index = index + N;
        pulse = data(start_index:end_index,1)';

        % MS Clutter Rejection
        pulse = pulse - mean(pulse);    % Subtract column mean to each column
        
        sfft_ms = 20*log10(abs(fft(pulse, 4*N)));  % dft using zero padding
        sfft_ms = sfft_ms(1:2*N) - max(sfft_ms); % Normalize data  
        
        % update spectrum
        spec_range(2:end,:) = spec_range(1:end-1,:);
        spec_range(1,:) = sfft_ms;
        set(spectrogram, 'CData', spec_range); pause(0.05)
            
        index = index + N;
        continue
    end

    index = index + 1;

%     while upchrip == true
%         if sync(index+1,1) > 0 && sync(index,1) < 0
%             end_index = index;
%             pulse = data(start_index:end_index,1);
% 
%             % fft and 
%         
%             index = index + 1;
%             
%             upchirp = false;
%         end
%     end
end

