close all; clear; clc;

%% COMPUTE VELOCITY FROM AUDIO FILE
% Read the audiofile 
% [y,Fs] = audioread('Velocity_Test_File.m4a'); 
[y,Fs] = audioread('audacity_recordings\multiple_targets_velocity.wav'); 


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
f2_smooth = smoothdata(f2,1);
% f2_smooth = f2;
vel_max1 = zeros(size(f2_smooth,1),1);
vel_max2 = zeros(size(f2_smooth,1),1);
for i = 1:size(f2_smooth,1)
    [~, locs, ~, prominence] = findpeaks(f2_smooth(i,:),'MinPeakHeight', -20);
    [~, ind] = maxk(prominence, 2);
    vel_max1(i) = vel2(locs(ind(1)));
    vel_max2(i) = vel2(locs(ind(2)));
end

% acc2 = zeros(size(vel_max2,1)-1,1);
% acc1 = zeros(size(vel_max1,1)-1,1);
% 
% 
% for i = 2:size(vel_max2,1)
%     der = abs((vel_max2(i) - vel_max2(i-1))/(time2(2)-time2(1)));
%     acc2(i-1) = der;
%     if acc2(i-1) > 10
%         vel_max2(i) = 0;
%     end
% end
% 
% for i = 2:size(vel_max1,1)
%     der = abs((vel_max1(i) - vel_max1(i-1))/(time1(2)-time1(1)));
%     acc1(i-1) = der;
%     if acc1(i-1) > 10
%         vel_max1(i) = 0;
%     end
% end
% 
% vel_max1(1) = 0;
% vel_max2(1) = 0;

figure(2);
%[~, f_ind] = maxk(f2,2,2);
%vel_max1 = vel2(f_ind(:,1));
%vel_max2 = vel2(f_ind(:,2));
hold on
plot(time2, vel_max1)
plot(time2, vel_max2)
% plot(time2(2:end), acc2)
% plot(time2(2:end), acc1)
grid on
% ylim([0, 7]);
title("Speed-time plot")
xlabel("time"), ylabel("Speed")
legend('vel1','vel2')
hold off

