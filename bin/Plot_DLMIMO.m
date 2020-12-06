%%
%% Plat_SUSISO.m - A MATLAB script code for plotting a C++ SUSISO result that is BLER and Throughput.
%%
%% Copyright(C) 2018 KAIST. All rights reserved. 
%%--------------------------------------------------------------------------------------------------------------------
CQI = 1:1;
SNR_stepsize_num = 4;   % 20;

%% SNR
SNR_start = [ 16.0, -14.0, -13.0, -10.0, -7.0, -5.0, -3.0, -1.0,  1.0,  4.0,  6.0,  8.0, 10.0, 12.0];
SNR_end   = [ 26.0,  10.0,  11.0,  12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0] + 30.0;

SNR_stepsize =  (SNR_end - SNR_start) / (SNR_stepsize_num - 1);
SNR_range = zeros(SNR_stepsize_num, length(SNR_start));
for ind_cqi = 1:length(SNR_start)
    SNR_range(:,ind_cqi) = SNR_start(ind_cqi) : SNR_stepsize(ind_cqi)  : SNR_end(ind_cqi);
end


BLER = dlmread('OUTPUT_CPP_BLER.csv');
Throughput = dlmread('OUTPUT_CPP_Throughput.csv');

legend_set = cell(length(CQI), 1);
for ind_cqi = 1:length(CQI)
    legend_set{ind_cqi} = ['CQI' num2str(CQI(ind_cqi))];
end

figure;
axes('YScale', 'log');
hold all; grid on; box on;
for ind_cqi = 1:length(CQI)
    plot(SNR_range(:,CQI(ind_cqi)), BLER(ind_cqi,:), 'LineWidth', 2);
end
legend(legend_set, 'Location', 'EastOutside');
xlabel('SNR [dB]');
ylabel('BLER');
title(' BLER DL_MIMO (CDL_A)');
ylim([0.01 1]);

figure;
hold all; grid on; box on;
for ind_cqi = 1:length(CQI)
    plot(SNR_range(:,CQI(ind_cqi)), Throughput(ind_cqi,:), 'LineWidth', 2);
end
legend(legend_set, 'Location', 'EastOutside');
xlabel('SNR [dB]');
ylabel('Throughput [Mbps]');
title(' Throughput DL_MIMO (CDL_A)');