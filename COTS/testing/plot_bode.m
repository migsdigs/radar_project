T = readtable('cots_bode.xlsx');
Tmat = T{:,:};

cutoff = 15.915; % [kHz]

% Oscilloscope
freq_o = Tmat(15,1);
ampl_o = Tmat(15,7);
str_o = " "+freq_o+"kHz";
plot(Tmat(2:end,1), Tmat(2:end,7)); grid on; ylabel("dB"); xlabel("kHz"); hold on;
plot(freq_o, ampl_o, '.', 'MarkerSize', 10); text(freq_o, ampl_o, str_o); hold on;

% Multimeter
freq_m = Tmat(16,13);
ampl_m = Tmat(16,19);
str_m = " "+freq_m+"kHz";
plot(Tmat(2:end,13), Tmat(2:end,19)); grid on; ylabel("dB"); xlabel("kHz"); hold on;
plot(freq_m, ampl_m, '.', 'MarkerSize', 10); text(freq_m, ampl_m, str_m); hold on;

legend('Oscilloscope','','Multimeter','') 