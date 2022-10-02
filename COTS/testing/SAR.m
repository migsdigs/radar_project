close all; clear; clc;

% Read the audiofile 
[y,Fs] = audioread('SAR_Test_File.m4a'); 

%% SETTING UP DATA MATRIX
% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = -y(:,2); % Sync data (square waveform)

Trp = 0.25;         % [s] Duration of measuring at each position
Nrp = Trp*Fs;       % [ ] Number of samples at each position
Tp = 20e-3;         % [s] Upchirp length
N = Tp*Fs;          % [ ] Number of samples per upchirp
c = 299792458;      % [m/s] Speed of light
fc = 2.43e9;        % [Hz] Carrier frequency
lambda = c/fc;      % [m] WaveLength
delta_x = lambda/2;
fStart = 2.408e9;   % Start Frequency [Hz]
fStop = 2.495e9;    % Stop Frequency  [Hz]


% sync = (sync > 0.25)'; % Set sync signal to 0 or 1
% 
% % Find start of positions. This is where there is at least N+1 zeros
% % followed by N ones
% sequenceToFind = [zeros(1,2*N), ones(1, 0.5*N)];
% found = strfind(sync, sequenceToFind);
% startIndices = found + 2*N + Nrp;
% numberOfPositions = length(startIndices);
% 
% dataMatrix = zeros(numberOfPositions,Nrp);
% syncMatrix = zeros(numberOfPositions,Nrp);
% for k = 1:numberOfPositions
%     dataMatrix(k,:) = data(startIndices(k):startIndices(k)+Nrp-1);
%     syncMatrix(k,:) = sync(startIndices(k):startIndices(k)+Nrp-1);
% end
% 
% % APPLYING HILBERT TRANSFORM
% % First, integrate over all up-chirps
% dataIntegrated = zeros(numberOfPositions,N);
% sequenceToFind = [0, 0, 1, 1];
% for k = 1:numberOfPositions
%     upchirpStart = strfind(syncMatrix(k,:),sequenceToFind);
%     for j = 1:(length(upchirpStart)-1)
%         dataIntegrated(k,:) = dataIntegrated(k,:) + dataMatrix(k,(upchirpStart(j)+2):(upchirpStart(j)+1+N));
%     end
%     dataIntegrated(k,:) = dataIntegrated(k,:)/length(upchirpStart);
% end

% separate square waveform and data
data_bs = data;
data_sc = sync >= 0.25;
data_N = round(Tp * Fs);
data_Nrp = ceil(Trp * Fs);

% find begginings of positions
data_pos_dif = diff([false; ~data_sc; false]);
data_pos_len = find(data_pos_dif<0) - find(data_pos_dif>0);
data_pos_idx = find(data_pos_dif>0);

data_pos_valid = data_pos_len >= 2*data_Nrp; %2.5*con_Tp*data_Fs;
data_pos_idx = data_pos_idx(data_pos_valid);
data_pos_len = data_pos_len(data_pos_valid);
data_pos_idx = data_pos_idx + data_pos_len;
clear data_pos_dif data_pos_valid data_pos_len;

% add the guard band
%data_pos_idx = data_pos_idx + data_Nrp + 0.5*data_N;
%figure(); plot(data_sc); hold on; xline(data_pos_idx, 'r');

% assemble 2D matrices
data_pos_idx = cell2mat(arrayfun(@(x) x:x+data_Nrp-1, data_pos_idx, 'UniformOutput', false));
data_bs = data_bs(data_pos_idx);
data_sc = data_sc(data_pos_idx);
clear data_pos_idx;

% set number of possitions
data_P = size(data_bs, 1);

% process every position separatelly
sar_data = NaN(data_P, data_N);
for i = 1:size(data_bs,1)
    % find upchrip begginings
    upchrip_idx = diff([false, data_sc(i, :), false]);
    upchrip_len = find(upchrip_idx<0) - find(upchrip_idx>0);
    upchrip_idx = find(upchrip_idx>0);
    upchrip_idx = upchrip_idx(upchrip_len >= 0.9 * data_N);

    % average upchirp
    upchrip = cell2mat(arrayfun(@(x) data_bs(i, x:x+data_N-1), upchrip_idx, 'UniformOutput', false));
    sar_data(i, :) = mean(reshape(upchrip', data_N, length(upchrip_idx))', 1); %#ok<UDIM> 
end
clear upchrip_idx upchrip data_bs;

dataIntegrated = sar_data;

%Now for Hilbert transform:
%FFT each row
dataFFT = fft(dataIntegrated,N,2);
%IFFT positive frequencies
dataFFT = dataFFT(:,1:end/2);
dataIFFT = ifft(dataFFT,N,2);
%set NaN values to 1e-30
dataIFFT(isnan(dataIFFT)) = 1e-30;

%% APPLYING HANN WINDOWING
numberOfPositions = size(dataIntegrated,1);
windowed = zeros(numberOfPositions,N);
hanningWindow = hann(N)';
for k = 1:numberOfPositions
    windowed(k,:) = dataIFFT(k,:).*hanningWindow;
end

kr = linspace(4*pi/c*fStart,4*pi/c*fStop,N);
xa = linspace(-(delta_x*numberOfPositions)/2,(delta_x*numberOfPositions)/2,numberOfPositions);
figure(), clf()
imagesc(kr,xa,angle(windowed))
colorbar;

%% RANGE MIGRATION ALGORITHM
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

d_range_1 = 1;
d_range_2 = 100;
c_range_1 = -25; 
c_range_2 = 25;  
flipped = fliplr(rot90(ifft_interp_dm)); 
d_index1 = round((size(flipped,1)/r_max) * d_range_1 * 4);
d_index2 = round((size(flipped,1)/r_max) * d_range_2 * 4);
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