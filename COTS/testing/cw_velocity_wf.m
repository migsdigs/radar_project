close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE


% Read the audiofile 
[y,Fs] = audioread('audacity_recordings\WavFile1.wav'); 

% Parameters and Constants
c = 299792458;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]

pad = 4;
wf = 64;

% Take care of data inversion by the sound card
data = -y(:,1); % Intensity of the received signal

% spec_vel = ones(wf, vel_sps) * -Inf;
% data_v = linspace(0, velo_max, vel_sps);
% data_t = linspace(0, Tp*wf,wf);

spec_vel = ones(wf, N) * -Inf;

delta_f = linspace(0, Fs/2, N); 
data_v = (delta_f * c)/(2 * f_center);
data_t = linspace(0, wf*Tp, wf);


% create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
spectrogram = imagesc(data_v, data_t, spec_vel, [-10 0]); colorbar; set(gca,'XLim',[0 20]);
title("Velocity vs Time Spectrogram"); xlabel('Velocity [m/s]'); ylabel('Time [s]');

for i = 1 : floor(size(data,1)/N)
    % MS clutter rejection
    pulse = data((N*(i-1)+1) : i*N,1)';
    pulse = pulse - mean(pulse);

    % fft and normalization
    sfft = 20*log10(abs(fft(pulse, N*4)));
%     sfft = sfft(1,1:size(sfft,2)/2);
    sfft = sfft(1,1:N);
    sfft = sfft - max(sfft);


    % update spectrum
    spec_vel(2:end,:) = spec_vel(1:end-1,:);
    spec_vel(1,:) = sfft;
    set(spectrogram, 'CData', spec_vel); pause(0.05)
    
end

