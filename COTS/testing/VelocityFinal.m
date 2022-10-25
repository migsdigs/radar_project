clear, clc

% Read the audiofile 
[y,Fs] = audioread('audacity_recordings\CW_MULTI_REDO.wav'); 
numRidges = 2; % Number of targets

% Take care of data inversion by the sound card
data = -y(:,1);

% In case of IQ Demodulation, choose which data to use from here
% [I,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_REAL.wav'); 
% [Q,Fs] = audioread('audacity_recordings/SDR_CWIF_BREATHING_MIGUEL_IMAG.wav'); 
% % Composite the complex value and its conjugate
% data = complex(I,Q);
% data = conj(data);

% Parameters
c = 299792458;                % Speed of light [m/s]
f_center = 2.45e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

fridgeLimit = 10;
fridgeLength = round((2*f_center*fridgeLimit*4)/c); % To limit the search of tfridge
freqs = linspace(fridgeLimit, 0, fridgeLength);


% Parse the data
X = mod(-mod(length(data), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';

% MS Clutter Rejection
final_data = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract column mean to each column

% FFT 
f = abs(fft(final_data, 4*N, 2));
f = 20*log10(f);
f = f(:,1:size(f, 2) / 2);

% Normalization 1
f1 = f - max(max(f));
delta_f1 = linspace(0, Fs/2, size(f1, 2)); 
vel1 = (delta_f1 * c)/(2 * f_center);
time1 = linspace(1, Tp * size(f1, 1), size(f1, 1));

% Normalization 2
f_row_max = max(f, [], 2);
f2 = f - f_row_max;
delta_f2 = linspace(0, Fs/2, size(f2, 2)); 
vel2 = (delta_f2 * c)/(2 * f_center);
time2 = linspace(1, Tp * size(f2,1), size(f2,1));

% find fridges
sfft_fridge = fft(final_data, 4*N, 2);
[fridge, ~, ~] = tfridge(rot90(sfft_fridge(:, 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',10);
% [fridge1, ~, ~] = tfridge(rot90(sfft_fridge(1:floor(end*7/16), 1:fridgeLength)), freqs, 0.3,'NumRidges',numRidges,'NumFrequencyBins',10);
% [fridge2, ~, ~] = tfridge(rot90(sfft_fridge((floor(end*7/16)+1):floor(end*9/16), 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',3);
% [fridge3, ~, ~] = tfridge(rot90(sfft_fridge((floor(end*9/16)+1):floor(end*13/16), 1:fridgeLength)), freqs, 0.2,'NumRidges',numRidges,'NumFrequencyBins',10);
% [fridge4, ~, ~] = tfridge(rot90(sfft_fridge((floor(end*13/16)+1):end, 1:fridgeLength)), freqs, 1,'NumRidges',numRidges,'NumFrequencyBins',1);
% fridge = [fridge1; fridge2; fridge3; fridge4];
% Convert to velocity
velExtracted = 160*(fridge * c) / (2*f_center);


%% Plotting
figure();
sgtitle('Multiple target CW velocity measurement');

spectogram=subplot(1,2,1);
imagesc(vel1, time1, f1);
caxis([-35 0]);
colorbar;
set(gca,'XLim',[0 10]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("Velocity spectogram");

subplot(1,2,2)
rangePlot = plot(time1,velExtracted);
title("Velocity vs Time plot"); xlabel('Time [s]'); ylabel('Velocity [m/s]');
