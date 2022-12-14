close all; clear; clc;

% Read the audiofile 
[y,Fs] = audioread('SAR_Test_File.m4a'); 

%% SETTING UP DATA MATRIX
% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

Trp = 0.25;         % Duration of measuring at each position [s]
Nrp = Trp*Fs;       % Number of samples at each position
Tp = 20e-3;         % Upchirp length [s]
N = Tp*Fs;          % Number of samples per upchirp
c = 299792458;      % Speed of light [m/s]
fc = 2.43e9;        % Carrier frequency [Hz]
lambda = c/fc;      % WaveLength [m]
delta_x = lambda/2;
fStart = 2.408e9;   % Start Frequency [Hz]
fStop = 2.495e9;    % Stop Frequency  [Hz]

% separate square waveform and data
data_backscatter = data;
data_sync = sync >= 0.25;
data_N = round(Tp * Fs);
data_Nrp = ceil(Trp * Fs);

% find begining of positions
data_pos_dif = diff([false; ~data_sync; false]);
data_pos_len = find(data_pos_dif<0) - find(data_pos_dif>0);
data_pos_idx = find(data_pos_dif>0);

data_pos_valid = data_pos_len >= 2*data_Nrp; %2.5*con_Tp*data_Fs;
data_pos_idx = data_pos_idx(data_pos_valid);
data_pos_len = data_pos_len(data_pos_valid);
data_pos_idx = data_pos_idx + data_pos_len-1;
clear data_pos_dif data_pos_valid data_pos_len;

data_pos_idx = cell2mat(arrayfun(@(x) x:x+data_Nrp-1, data_pos_idx(1:end-1), 'UniformOutput', false));
data_backscatter = data_backscatter(data_pos_idx);
data_sync = data_sync(data_pos_idx);
clear data_pos_idx;

% number of positions
data_P = size(data_backscatter, 1);

% positions processed
sar_data = NaN(data_P, data_N);
for i = 1:size(data_backscatter,1)
    % find upchrip begginings
    upchrip_idx = diff([false, data_sync(i, :), false]);
    upchrip_len = find(upchrip_idx<0) - find(upchrip_idx>0);
    upchrip_idx = find(upchrip_idx>0);
    upchrip_idx = upchrip_idx(upchrip_len >= 0.9 * data_N);

    % integrate upchirps
    upchrip = cell2mat(arrayfun(@(x) data_backscatter(i, x:x+data_N-1), upchrip_idx, 'UniformOutput', false));
    sar_data(i, :) = mean(reshape(upchrip', data_N, length(upchrip_idx))', 1); %#ok<UDIM> 
end
clear upchrip_idx upchrip data_backscatter;

dataIntegrated = sar_data;

% HILBERT TRANSFORM
% FFT each row
dataFFT = fft(dataIntegrated,N,2);

% IFFT positive frequencies
dataFFT = dataFFT(:,1:end/2);
dataIFFT = ifft(dataFFT,N,2);

% set NaN=1e-30
dataIFFT(isnan(dataIFFT)) = 1e-30;

% APPLYING HANN WINDOWING
numberOfPositions = size(dataIntegrated,1);
windowed = zeros(numberOfPositions,N);
hanningWindow = hann(N)';
for k = 1:numberOfPositions
    windowed(k,:) = dataIFFT(k,:).*hanningWindow;
end
% windowed = dataIFFT;

kr = linspace(4*pi/c*fStart,4*pi/c*fStop,N);
xa = linspace(-(delta_x*numberOfPositions)/2,(delta_x*numberOfPositions)/2,numberOfPositions);
figure(), clf()
imagesc(kr,xa,angle(windowed))
colorbar;

% RANGE MIGRATION ALGORITHM
% Add zero-padding
padding = 2048;
temp_zeros = zeros(padding, size(windowed,2));
for i = 1:size(windowed,2)
    index = round((padding-size(windowed,1))/2);
    temp_zeros(index+1:(index + size(windowed,1)),i) = windowed(:,i);
end
windowed_padded = temp_zeros;

data_matrix = fftshift(fft(windowed_padded,[],1),1);
kx = linspace((-2*pi/lambda),(2*pi/lambda),(size(data_matrix,1)));

% Plot the magnitude in dB of cross-range FFT
figure();
imagesc(kr, kx, 20*log10(abs(data_matrix)), [max(max(20*log10(abs(data_matrix))))-40,max(max(20*log10(abs(data_matrix))))]);
xlabel('K_r (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

% Plot the phase of cross-range FFT
figure();
imagesc(kr, kx, angle(data_matrix));
xlabel('K_r (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

% Stolt interpolation
ky = zeros(padding, N);
interpolated_dm = zeros(padding, padding/2);
for i = 1:padding
    ky(i,:) = sqrt(kr.^2 - kx(i)^2);
end

k_yel = floor(min(min(ky)));
k_yezpad = ceil(max(max(ky)));
k_ye = linspace(k_yel, k_yezpad, padding/2);

for i = 1:padding
    interpolated_dm(i,:) = (interp1(ky(i,:), data_matrix(i,:), k_ye));
end

interpolated_dm(isnan(interpolated_dm)) = 1e-30;

% Plot the phase after Stolt interpolation
figure();
imagesc(k_ye, kx, angle(interpolated_dm)); 
xlabel('K_y (rad/m)'); ylabel('K_x (rad/m)');
colorbar;

ifft_interp_dm = ifft2(interpolated_dm, (padding*4), ((padding/2)*4));

dfy = c*(k_yezpad - k_yel)/(4*pi);
r_max = ((4*c*padding/2)/(2*dfy));
rail_r_max = padding * delta_x;

constant = 3;

d_range_1 = 1;
d_range_2 = 100;
c_range_1 = -25; 
c_range_2 = 25;  
flipped = fliplr(rot90(ifft_interp_dm)); 
d_index1 = round((size(flipped,1)/r_max) * d_range_1 * constant);
d_index2 = round((size(flipped,1)/r_max) * d_range_2 * constant);
c_index1 = round((size(flipped,2)/rail_r_max) * (c_range_1+(rail_r_max/2)));
c_index2 = round((size(flipped,2)/rail_r_max) * (c_range_2+(rail_r_max/2)));

truncated = flipped(d_index1:d_index2, c_index1:c_index2);
down_range = linspace(-d_range_1, -d_range_2, size(truncated,1));
cross_range = linspace(c_range_1, c_range_2, size(truncated, 2));

truncated_final = truncated;
for j = 1:size(truncated,2)
    truncated_final(:,j) = (truncated(:,j)').*(abs(down_range)).^2;
end
truncated_db = 20 * log10(abs(truncated_final));

figure();
imagesc(cross_range, down_range, truncated_db);  
xlim([-70 70]); caxis([-60 -15]); colorbar;