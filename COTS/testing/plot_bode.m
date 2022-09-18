T = readtable('cots_bode.xlsx');
Tmat = T{:,:};

cutoff = 15.915; % [kHz]

% Oscilloscope
freq_o = Tmat(15,1);
ampl_o = Tmat(15,7);
y_lim = -5;
x_lim = 2;
str_o = " "+freq_o+"kHz";
plot(Tmat(2:end,1), Tmat(2:end,7), 'Color', 'blue'); 
hold on;
plot(freq_o, ampl_o, '.', 'MarkerSize', 10); text(freq_o, ampl_o, str_o); 
hold on;
line([freq_o freq_o]', [y_lim ampl_o]', 'Color', 'blue', 'LineStyle','--');
line([x_lim freq_o]', [ampl_o ampl_o]', 'Color', 'blue', 'LineStyle','--');
hold on;

% Multimeter
freq_m = Tmat(16,13);
ampl_m = Tmat(16,19);
str_m = " "+freq_m+"kHz";
plot(Tmat(2:end,13), Tmat(2:end,19), 'Color', 'red');
hold on;
plot(freq_m, ampl_m, '.', 'MarkerSize', 10); text(freq_m, ampl_m, str_m); 
hold on;
line([freq_m freq_m]', [y_lim ampl_m]', 'Color', 'red', 'LineStyle','--');
line([x_lim freq_m]', [ampl_m ampl_m]', 'Color', 'red', 'LineStyle','--');

title("3dB roll-offs"); grid on; ylabel("dB"); xlabel("kHz"); ylim([y_lim 1]); xlim([x_lim 20]); 
legend('Oscilloscope','','','','Multimeter','','','') 