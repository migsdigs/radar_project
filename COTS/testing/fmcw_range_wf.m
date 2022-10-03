close all; clear; clc;
%% NEW SCRIPT
% Read the audiofile 
[y,Fs] = audioread('audacity_recordings\single_target_range.wav'); 

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

wf = 512;                     % Waterfall steps
r_max = 20;                   % Max range
v_max = 5;                    % Max velocity
padding = 4;

range_sps = round((2*bandwidth*r_max*padding)/c);
velocity_sps = round((2*f_center*v_max)/c);

% prepare spectrum, sweep history and chirp spectra arrays
spec = ones(wf, range_sps) * -Inf;
sweep = zeros(3, 2*N);

sfft_up = complex(zeros(wf, max(range_sps,velocity_sps)), 0);
sfft_dn = complex(zeros(wf, max(range_sps,velocity_sps)), 0);

% compute time and range vectors
data_d = linspace(0, r_max, range_sps);
data_t = linspace(0, 2 * Tp * wf, wf);

chrip_f = linspace(Fs / (padding / 2), 0, padding*N);
chrip_f = chrip_f(1:max(range_sps,velocity_sps));
chrip_t = linspace(0, 2 * Tp * wf, wf);

% create figures
figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
subplot(2,1,1)
velocityPlot = plot(data_t,NaN(1,wf));
title("Velocity vs Time plot"); xlabel('Time [s]'); ylabel('Velocity [m/s]');

subplot(2,1,2)
rangePlot = plot(data_t,NaN(1,wf));
title("Range vs Time plot"); xlabel('Time [s]'); ylabel('Range [m]');

i = 1;

while i < (size(data, 1) - 4*N)
    % update sweep history and spectrum
    sweep(2:3, :) = sweep(1:2, :);
    spec(2:end, :) = spec(1:end-1, :);

    % update sfft history
    sfft_up(2:end, :) = sfft_up(1:end-1, :);
    sfft_dn(2:end, :) = sfft_dn(1:end-1, :);

    % separate square waveform and data
    i = i + find(sync(i:end) >= 0.25, 1, 'first');
    sweep(1, 1:N) = data(i:i+N-1);
    i = i + find(sync(i:end) <= -0.25, 1, 'first');
    sweep(1, N+1:end) = data(i:i+N-1);
    i = i + N;

    % 2-pulse MTI
    sweep_mti = sweep(1, :) - sweep(2, :); % - (sweep(2, :) - sweep(3, :));

    % separate FFFT for each churp
    sfft_tmp = ifft(sweep_mti(1:data_N), padding*N);
    sfft_up(1,:) = sfft_tmp(1:max(range_sps,velocity_sps));
    sfft_tmp = ifft(sweep_mti(N+1:end), padding*N);
    sfft_dn(1,:) = sfft_tmp(1:max(range_sps,velocity_sps));

    % find fridges
    [fridge_up, ~, ~] = tfridge(rot90(sfft_up(:, 1:max(range_sps,velocity_sps))), chrip_f, 1);
    [fridge_dn, ~, ~] = tfridge(rot90(sfft_dn(:, 1:max(range_sps,velocity_sps))), chrip_f, 1);
 
    % compute range and velocity
    f_beat = (fridge_dn + fridge_up) / 2;
    rng = (f_beat * c * Tp) / bandwidth;
    rng = (rng - 1485);
    f_dopler = (fridge_dn - fridge_up) / 2;
    vel = (c * f_dopler) / (2 * f_center);
    vel = vel - 0.2;
    
    % update spectrum and all plots
    spec(1,:) = 20*log10(abs(sfft_up(1, 1:range_sps)));
    spec(1,:) = spec(1,:) - max(spec(1,:));
    
    %set(sp, 'CData', spec);
    set(rangePlot, 'YData', rng);
    set(velocityPlot, 'YData', vel);
    
    pause(0.01);
end


%% COMPUTE RANGE FROM AUDIO FILE
% Read the audiofile 
[y,Fs] = audioread('audacity_recordings\single_target_range.wav'); 

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
dr = c / (2*bandwidth);       % Range resolution [m]
Tp = 20e-3;                   % Pulse width (up-chirp) [s]
N = Tp * Fs;                  % Number of samples per up-chrip pulse
max_range = (N * dr)/2;       % Maximum range [m]

wf = 512;
usingDoppler = false;

% % plotting
spec_range = ones(wf, 2*N) * -Inf;
spec_velocity = ones(wf, 2*N) * -Inf;

range = linspace(0, max_range, N*2); 
delta_f1 = linspace(0, Fs/2, 2*N); 
vel1 = (delta_f1 * c)/(2 * f_center);
velocity = zeros(wf,1);
rangePeaks = zeros(wf,1);
data_t = linspace(0, wf*Tp*2, wf);

% create figure
figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
spectrogramRange = imagesc(range, data_t, spec_range, [-20 0]); colorbar; set(gca,'XLim',[0 100]);
title("Range vs Time Spectrogram"); xlabel('Range [m]'); ylabel('Time [s]');

figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
spectrogramVelocity = imagesc(vel1, data_t, spec_velocity, [-20 0]); colorbar; set(gca,'XLim',[0 100]);
title("Velocity vs Time Spectrogram"); xlabel('Velocity [m/s]'); ylabel('Time [s]');

figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
velocityPlot = plot(data_t,velocity);
title("Velocity vs Time plot"); xlabel('Time [s]'); ylabel('Velocity [m/s]');

figure('Units', 'normalized', 'OuterPosition', [0 0 0.4 0.4]);
rangePlot = plot(data_t,rangePeaks);
title("Range vs Time plot"); xlabel('Time [s]'); ylabel('Range [m]');

index = 1;
start_index = 1;
end_index = 0;
upchrip = false;
pulse = zeros(1,[]);
oldpulse = 0;
x = true;

while x == true
    if sync(index+1,1) > 0 && sync(index,1) < 0
        if usingDoppler
            start_index = index + 1;
            end_index = index + 2*N;
            pulse = data(start_index:end_index,1)';
            
            newpulse = pulse - oldpulse; %2-point MTI
            oldpulse = pulse;
            pulse = newpulse;

            [rangeNew, velocityNew] = FMCW(pulse);
%             velocity(2:end) = velocity(1:end-1);
%             velocity(1) = velocityNew;
%             rangePeaks(2:end) = rangePeaks(1:end-1);
%             rangePeaks(1) = rangeNew;
%             
%             set(velocityPlot, 'YData', velocity);
%             set(rangePlot, 'YData', rangePeaks);
            spec_range(2:end,:) = spec_range(1:end-1,:);
            spec_range(1,:) = rangeNew;
            set(spectrogramRange, 'CData', spec_range);

            spec_velocity(2:end,:) = spec_velocity(1:end-1,:);
            spec_velocity(1,:) = velocityNew;
            set(spectrogramVelocity, 'CData', spec_velocity);

        else
            start_index = index + 1;
            end_index = index + N;
            pulse = data(start_index:end_index,1)';
            
            newpulse = pulse - oldpulse; %2-point MTI
            oldpulse = pulse;
            pulse = newpulse;

            % MS Clutter Rejection
            pulse = pulse - mean(pulse);    % Subtract column mean to each column
        
            sfft_ms = 20*log10(abs(fft(pulse, 4*N)));  % dft using zero padding
            sfft_ms = sfft_ms(1:2*N) - max(sfft_ms); % Normalize data  

            % update range spectrum
            spec_range(2:end,:) = spec_range(1:end-1,:);
            spec_range(1,:) = sfft_ms;
            set(spectrogramRange, 'CData', spec_range);

            %update velocity specturm
            velocity(2:end) = velocity(1:end-1);
            [~, ind1] = max(sum(spec_range(1:10,:),1));
            [~, ind2] = max(sum(spec_range(20:30,:),1));
            velocity(1) = (ind1 - ind2)*max_range/(20*2*Tp*2*N);
            set(velocityPlot, 'YData', velocity);
        end
        pause(0.01);    
        index = index + N;
        continue
    end

    index = index + 1;

end


function [range, velocity] = FMCW(data)
    Fs = 44100;
    f_center = 2.43e9;            % Center Frequency [Hz]
    c = 299792458;                % Speed of light [m/s]
    f_start = 2.408e9;            % Start Frequency [Hz]
    f_stop = 2.495e9;             % Stop Frequency  [Hz]
    BW = f_stop - f_start;        % [Hz]
    Tp = 20e-3;                   % Pulse width [s]
    N = Tp * Fs;                  % Number of samples per pulse

    % Parse up-chirp data according to the sync data
%     sync_pulse = zeros(length(sync),1);
%     sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1

    % Data Upchirp
%     for i = 2:(size(sync_pulse)-N)
%         if sync_pulse(i) ~= sync_pulse(i-1)
%             %There is a 0.25% chance that the change is exactly at the first element,
%             %then it will give an error
%             if sync_pulse(i) == 1
%                 up_data_parsed = data(i:(i+N-1));
%                 down_data_parsed = [data(1:i);data((i+N-1):end)];
%             else
%                 up_data_parsed = [data(1:i);data((i+N-1):end)];
%                 down_data_parsed = data(i:(i+N-1));
%             end
%         end
%     end
    up_data_parsed = data(1:end/2);
    down_data_parsed = data((end/2+1):end);
    
    % Remove noise
    up_data = up_data_parsed - mean(up_data_parsed);
    down_data = down_data_parsed - mean(down_data_parsed);

    % Do FFT
    fft_up = fft(up_data,4*N);
    fft_down = fft(down_data,4*N);

    % Get the beat frequencies
%     [~, ind_up] = max(fft_up);
%     [~, ind_down] = max(fft_down);
%     freq_up = 1/(2*Tp)*ind_up/length(fft_up);
%     freq_down = 1/(2*Tp)*ind_down/length(fft_down);
% 
%     range_freq = (freq_up + freq_down)/2;
%     velocity_freq = (freq_up-freq_down)/2;
% 
%     range = Tp*c*range_freq/(2*BW);
%     velocity = (velocity_freq*c)/(2*f_center);

    
    range = (fft_up + fft_down)/2;
    range = 20*log10(abs(range));
    range = range(1:2*N) - max(range); %Normalise data

    velocity = (fft_up - fft_down)/2;
    velocity = 20*log10(abs(velocity));
    velocity = velocity(1:2*N) - max(velocity); %Normalise data
end
