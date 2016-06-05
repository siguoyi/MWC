clear all;
close all;

ber = xlsread('ber.xlsx');
bb = ones(9,1);
meanBer = mean(ber,2);
meanBer = bb - meanBer/64;
semilogy(meanBer, '-o','linewidth',2)
set(gca,'xticklabel',get(gca,'xTick')*5-15);
% grid on
xlabel('Sampling SNR');
ylabel('BER');

