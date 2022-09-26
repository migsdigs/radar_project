clear, clc, close all

recordingTime = 10;
saveAudio = true;
audiofileSaveName = 'RealTimeTest.wav';

rangeFlag = true; %false for velocity

if rangeFlag
    Tp = 20e-3; %Up-chirp time of FMCW
else
    Tp = 0.1; %For CW
end

% Parameters
Fs = 44100;
c = 299792458;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
N = Tp * Fs;                  % Number of samples per pulse
f_start = 2.408e9;            % Start Frequency [Hz]
f_stop = 2.495e9;             % Stop Frequency  [Hz]
bandwidth = f_stop - f_start; % [Hz]
dr = c / (2*bandwidth);       % Range resolution [m]
max_range = (N * dr)/2;       % Maximum range [m]

% Sound card parameters
nBits = "16-bit integer";
nChannels = 2;

% Audio reader system
deviceReader = audioDeviceReader(Fs,Fs*2*2*Tp,"BitDepth",nBits,"NumChannels",nChannels);
setup(deviceReader)
if saveAudio
    fileWriter = dsp.AudioFileWriter(audiofileSaveName,'FileFormat','WAV');
end

%% CW SETUP
if ~rangeFlag
    pad = 4;
    wf = 64;

    spec_vel = ones(wf, N) * -Inf;

    delta_f = linspace(0, Fs/2, N); 
    data_v = (delta_f * c)/(2 * f_center);
    data_t = linspace(0, wf*Tp, wf);

    % create figure
    figure('Units', 'normalized', 'OuterPosition', [0 0 0.8 0.8]);
    spectrogram = imagesc(data_v, data_t, spec_vel, [-10 0]); colorbar; set(gca,'XLim',[0 20]);
    title("Velocity vs Time Spectrogram"); xlabel('Velocity [m/s]'); ylabel('Time [s]');

end

%% FMCW SETUP
if rangeFlag
    wf = 128;

    spec_range = ones(wf, 2*N) * -Inf;

    range = linspace(0, max_range, N*2);
    data_t = linspace(0, wf*Tp*2*2, wf);

    % create figure
    figure('Units', 'normalized', 'OuterPosition', [0 0 0.8 0.8]);
    spectrogram = imagesc(range, data_t, spec_range, [-10 0]); colorbar; set(gca,'XLim',[0 30]);
    title("Range vs Time Spectrogram"); xlabel('Velocity [m/s]'); ylabel('Time [s]');

    index = 1;
    start_index = 1;
    end_index = 0;
    upchrip = false;
    pulse = zeros(1,[]);
    x = true;
    oldpulse = 0;
end


%% THIS IS KINDA FAILING I THINK
% To stop the recording, press a key
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
%set(a, 'KeyPressFcn', @myKeyPressFcn)

%% Starting the RT Radar
time = 0;
while ~KEY_IS_PRESSED
    acquiredAudio = deviceReader();
    if saveAudio
        fileWriter(acquiredAudio);
    end
    
    % Processing goes here
    data = -acquiredAudio(:,1); % Radar backscatter data (received reflected signal)
    sync = acquiredAudio(:,2); % Sync data (square waveform)  

    if rangeFlag
        found = find(sync(1:(end-1))<0 & sync(2:end) > 0);
        start_index = min(found);
            end_index = start_index + N - 1;
            pulse = data(start_index:end_index,1)';
            newpulse = pulse - oldpulse; %2-point MTI
            oldpulse = pulse;
            pulse = newpulse;
            
            % MS Clutter Rejection
            pulse = pulse - mean(pulse);    % Subtract column mean to each column
        
            sfft_ms = 20*log10(abs(fft(pulse, 4*N)));  % dft using zero padding
            sfft_ms = sfft_ms(1:2*N) - max(sfft_ms); % Normalize data  
            
            % update spectrum
            spec_range(2:end,:) = spec_range(1:end-1,:);
            spec_range(1,:) = sfft_ms;
            set(spectrogram, 'CData', spec_range);

    else
        pulse = data(1:N)';
        % MS clutter rejection
        pulse = pulse - mean(pulse);

        % fft and normalization
        sfft = 20*log10(abs(fft(pulse, N*4)));
        sfft = sfft(1,1:N);
        sfft = sfft - max(3, max(sfft));

        % update spectrum
        spec_vel(2:end,:) = spec_vel(1:end-1,:);
        spec_vel(1,:) = sfft;
        set(spectrogram, 'CData', spec_vel);
    end
    
    if ~isempty(get(figure(1),'CurrentCharacter'))
        break
    end
end

disp('Recording finished')
release(deviceReader)
if saveAudio
    release(fileWriter)
end

%% Functions to use
function [range, velocity] = FMCW(data, sync)
    Fs = 44100;
    f_center = 2.43e9;            % Center Frequency [Hz]
    c = 299792458;                % Speed of light [m/s]
    f_start = 2.408e9;            % Start Frequency [Hz]
    f_stop = 2.495e9;             % Stop Frequency  [Hz]
    BW = f_stop - f_start;        % [Hz]
    Tp = 20e-3;                   % Pulse width [s]
    N = Tp * Fs;                  % Number of samples per pulse

    % Parse up-chirp data according to the sync data
    sync_pulse = zeros(length(sync),1);
    sync_pulse(:,1) = (sync > 0);                  % Set sync square waveform between 0 and 1

    % Data Upchirp
    for i = 2:(size(sync_pulse)-N)
        if sync_pulse(i) ~= sync_pulse(i-1)
            %There is a 0.25% chance that the change is exactly at the first element,
            %then it will give an error
            if sync_pulse(i) == 1
                up_data_parsed = data(i:(i+N-1));
                down_data_parsed = [data(1:i);data((i+N-1):end)];
            else
                up_data_parsed = [data(1:i);data((i+N-1):end)];
                down_data_parsed = data(i:(i+N-1));
            end
        end
    end
    
    % Remove noise
    up_data = up_data_parsed - mean(up_data_parsed);
    down_data = down_data_parsed - mean(down_data_parsed);

    % Do FFT
    fft_up = fft(up_data,4*N);
    fft_down = fft(down_data,4*N);

    % Get the beat frequencies
    [~, ind_up] = max(fft_up);
    [~, ind_down] = max(fft_down);
    freq_up = 1/(2*Tp)*ind_up/length(fft_up);
    freq_down = 1/(2*Tp)*ind_down/length(fft_down);

    range_freq = (freq_up + freq_down)/2;
    velocity_freq = (freq_up-freq_down)/2;

    range = Tp*c*range_freq/(2*BW);
    velocity = (velocity_freq*c)/(2*f_center);    
end

function myKeyPressFcn(hObject, event)
    global KEY_IS_PRESSED
    KEY_IS_PRESSED  = 1;
end
