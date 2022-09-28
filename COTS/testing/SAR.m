close all; clear; clc;

% Read the audiofile 
[y,Fs] = audioread('SAR_Test_File.m4a'); 

%% SETTING UP DATA MATRIX
% Separate the sync data and radar backscatter data and take care of data
% inversion by the sound card
data = -y(:,1); % Radar backscatter data (received reflected signal)
sync = y(:,2); % Sync data (square waveform)

Trp = 0.25;      %[s] Duration of measuring at each position
Nrp = Trp*Fs;   %[ ] Number of samples at each position
Tp = 20e-3;     %[s] Upchirp length
N = Tp*Fs;      %[ ] Number of samples per upchirp
c = 299792458;  %[m/s] Speed of light
fc = 2.43e9;    %[Hz] Carrier frequency
lambda = c/fc;  %[m] WaveLength
fStart = 2.408e9;            % Start Frequency [Hz]
fStop = 2.495e9;             % Stop Frequency  [Hz]


sync = (sync > 0.25)'; % Set sync signal to 0 or 1

% Find start of positions. This is where there is at least N+1 zeros
% followed by N ones
sequenceToFind = [zeros(1,2*N), ones(1, 0.5*N)];
found = strfind(sync, sequenceToFind);
startIndices = found + 2*N + Nrp;
numberOfPositions = length(startIndices);

dataMatrix = zeros(numberOfPositions,Nrp);
syncMatrix = zeros(numberOfPositions,Nrp);
for k = 1:numberOfPositions
    dataMatrix(k,:) = data(startIndices(k):startIndices(k)+Nrp-1);
    syncMatrix(k,:) = sync(startIndices(k):startIndices(k)+Nrp-1);
end

%% APPLYING HILBERT TRANSFORM
% First, integrate over all up-chirps
dataIntegrated = zeros(numberOfPositions,N);
sequenceToFind = [0, 0, 1, 1];
for k = 1:numberOfPositions
    upchirpStart = strfind(syncMatrix(k,:),sequenceToFind);
    for j = 1:(length(upchirpStart)-1)
        dataIntegrated(k,:) = dataIntegrated(k,:) + dataMatrix(k,(upchirpStart(j)+2):(upchirpStart(j)+1+N));
    end
    dataIntegrated(k,:) = dataIntegrated(k,:)/length(upchirpStart);
end

%Now for Hilbert transform:
%FFT each row
dataFFT = fft(dataIntegrated,N,2);
%IFFT positive frequencies
dataFFT = dataFFT(:,1:end/2);
dataIFFT = ifft(dataFFT,N,2);
%set NaN values to 1e-30
dataIFFT(isnan(dataIFFT)) = 1e-30;

%% APPLYING HANN WINDOWING
windowed = zeros(numberOfPositions,N);
hanningWindow = hann(N)';
for k = 1:numberOfPositions
    windowed(k,:) = dataIFFT(k,:).*hanningWindow;
end

Kr = linspace(4*pi/c*fStart,4*pi/c*fStop,N);
Xa = linspace(-(lambda/2*numberOfPositions)/2,(lambda/2*numberOfPositions)/2,numberOfPositions);
figure(), clf()
imagesc(Kr,Xa,angle(windowed))
colorbar;
%% RANGE MIGRATION ALGORITHM
