close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
% [y,Fs] = audioread('Velocity_Test_File.m4a'); 
[y,Fs] = audioread('audacity_recordings/SDR_CWIF_MULTIPLE.wav'); 


% Take care of data inversion by the sound card
data = -y(:,1); % Intensity of the received signal

% Parameters
c = 299792458;                % Speed of light [m/s]
f_center = 2.43e9;            % Center Frequency [Hz]
Tp = 0.1;                     % Pulse width [s]
N = Tp * Fs;                  % Number of samples per pulse

% Parse the data
X = mod(-mod(length(data), N), N);      % Used to find the previous divisible value with respect to length(up_data)
data_cut = data((N-X+1):end);           % Remove the first elements so that we can reshape up_data
data_parsed = reshape(data_cut, N, [])';

% MS Clutter Rejection
final_data = bsxfun(@minus, data_parsed, mean(data_parsed, 2)); % Subtract column mean to each column
% final_data = data_parsed;

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

% Plotting
f_c_plot = f_center/1e9;
figure(1);
a1=subplot(1,2,1);
imagesc(vel1, time1, f1);
caxis([-40 0]);
colorbar;
set(gca,'XLim',[0 20]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("Normalization 1, Pulse time T_p="+Tp+"s, Center Frequency fc="+f_c_plot+"GHz");

a2=subplot(1,2,2);
imagesc(vel2, time2, f2);
caxis([-10 0]);
colorbar;
set(gca,'XLim',[0 20]);
xlabel('Velocity [m/sec]'); ylabel('Time [sec]');
title("Normalization 2, Pulse time T_p="+Tp+"s, Center Frequency fc="+f_c_plot+"GHz");


%% Velocity time plot
% First, find all peaks. Then check the peaks with biggest prominence and
% get the indices of these peaks in f2
k = 1; %Number of targets
treshhold = -15; %dB treshhold to indicate there is a target

f2_smooth = smoothdata(f2,1);
vel_max = zeros(size(f2_smooth,1),k);
for i = 1:size(f2_smooth,1)
    [~, locs, ~, prominence] = findpeaks(f2_smooth(i,:));
    [~, ind] = maxk(prominence, k);
    
    for j = 1:k
        % Only accept peaks if their intensity is above a treshhold
        if f2_smooth(i,locs(ind(j))) < treshhold
            vel_max(i,j) = 0;
        else
            vel_max(i,j) = vel2(locs(ind(j)));
        end
    end
end

figure(2), clf(), hold on
plot(time2, vel_max)
grid on
ylim([0, 7]);
title("Speed-time plot")
xlabel("time"), ylabel("Speed")
legend('vel1','vel2')
hold off

figure(3), clf();
subplot(1,3,1); plot(data); title("Data Captured");

subplot(1,3,2); imagesc(vel2, time2, f2); caxis([-10 0]); colorbar; set(gca,'XLim',[0 20]); 
xlabel('Velocity [m/s]'); ylabel('Time [s]'); title("Velocity vs Time Spectrogram"); 

subplot(1,3,3); plot(time2, vel_max); grid on; title("Velocity vs Time Plot");
ylim([0, 7]); xlabel("Time [s]"), ylabel("Velocity [m/s]"); legend('vel1','vel2')

