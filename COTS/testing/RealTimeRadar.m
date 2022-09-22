clear, clc, close all

recordingTime = 10;
saveAudio = true;
audiofileSaveName = 'RealTimeTest.wav';

% Parameters
Tp = 20e-3*2; %Up-chirp is 20ms

% Sound card parameters
Fs = 44100;
nBits = "16-bit integer";
nChannels = 2;

% Audio reader system
deviceReader = audioDeviceReader(Fs,Fs*Tp,"BitDepth",nBits,"NumChannels",nChannels);
setup(deviceReader)
if saveAudio
    fileWriter = dsp.AudioFileWriter(audiofileSaveName,'FileFormat','WAV');
end

velocityPlot = animatedline;
rangePlot = animatedline;

% To stop the recording, press a key
global KEY_IS_PRESSED
KEY_IS_PRESSED = 0;
gcf
set(gcf, 'KeyPressFcn', @myKeyPressFcn)

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
    
    %[range, velocity] = FMCW(data, sync);
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

    % Plotting goes here
    %addpoints(rangePlot, time, range)
    addpoints(velocityPlot, time, velocity);
    drawnow limitrate
    
    time = time + Tp;
    if ~isempty(get(figure(1),'CurrentCharacter'))
        break
    end
end
drawnow

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
