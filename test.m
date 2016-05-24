clear all;
close all;

tt = xlsread('tt.xlsx');
t = ones(1,20);
tt(1,:) = t - tt(1,:)/64;
figure(1)
subplot(211)
plot(tt(1,:), '-o')
title('BER')
subplot(212)
plot(tt(2,:), '-*b')
title('SNR')