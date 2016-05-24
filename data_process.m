clear all;
close all;

ber = xlsread('ber.xlsx');
bb = ones(9,1);
meanBer = mean(ber,2);
meanBer = bb - meanBer/64;
plot(meanBer, '-o','linewidth',2)
set(gca,'xticklabel',get(gca,'xTick')*5-15);
xlabel('Input SNR');
ylabel('BER');

