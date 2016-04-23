clear all;
close all;
clc;

%% Parameters

PPM = 20;
R = 1e5; % [bits/sec]���ı��ĸ�����
duration = 0.128; % [sec]
DataL = R*duration/2;  % Data length in symbols  �������ݳ���
NFFT = 128;  %��ʼֵ128

M = 4; % 4 QAM
snr = 20; % [dB]


Nsym = 4;           % Filter order in symbol durations  ��ʼֵΪ4
beta = 0.5;         % Roll-off factor  ����ϵ��
sampsPerSym = 10;    % Upsampling factor  ��Ƶϵ�����ɸı�Ƶ����С) ��ʼֵΪ10
L = sampsPerSym*Nsym + 1; % Raised cosine filter order  �������˲���
Fs = R * sampsPerSym;   % Sampling frequency ȡ��Ƶ��

cyclic_prefix_signal = [];

%% Raised cosine filter design

shape = 'Raised Cosine';  
% Specifications of the raised cosine filter with given order in symbols
rcosSpec = fdesign.pulseshaping(sampsPerSym, shape, 'Nsym,beta', Nsym, beta);
rcosFlt = design(rcosSpec);
rcosFlt.Numerator = rcosFlt.Numerator / max(rcosFlt.Numerator);

signal=(zeros(((NFFT+L)/NFFT)*DataL*sampsPerSym,1))';

fhz = [Fs/5,Fs/2.5];
for i=1:2

%% Tranceiver

bernoulli_binary_generator = randint(R*duration,1);

% Transmiter

bernoulli_two_samples = reshape(bernoulli_binary_generator, length(bernoulli_binary_generator)/2, 2);
dec = bi2de(bernoulli_two_samples,'left-msb'); % Bit to integer
modulated_data = qammod(dec,M); % 4 QAM  �����ϵ�ͼ

% figure(3)
% I = real(modulated_data);
% Q = imag(modulated_data);
% subplot(211);
% x=1:100;
% plot(x,Q(1:100),'rs');
% title('4QAM���ƺ���鲿');
% subplot(212);
% x=1:100;
% plot(x,I(1:100),'bd');
% title('4QAM���ƺ��ʵ��');

% Serial to Parallel ����������

cyclic_prefix_signal = [];
for id = 1:length(modulated_data)/NFFT
    ifft_signal = ifft(modulated_data((id-1)*NFFT+1:id*NFFT)); % ifft
    % adding cyclic prefix ����ѭ��ǰ׺
    cyclic_prefix = zeros(NFFT + L, 1);
    cyclic_prefix(1:L) = ifft_signal(end-L+1:end);
    cyclic_prefix(L+1:end) = ifft_signal;
    cyclic_prefix_signal = [cyclic_prefix_signal; cyclic_prefix];
end;  
%����ѭ��ǰ׺�Ժ���ź�
% figure(5)
% x=1:length(cyclic_prefix_signal);
% y1=abs(fft(cyclic_prefix_signal,length(cyclic_prefix_signal)));
% plot(x(1:length(cyclic_prefix_signal)/2),y1(1:length(cyclic_prefix_signal)/2));
% title('ifft�����ź�Ƶ��');

% D/A
signal_complex = filter(rcosFlt, upsample([cyclic_prefix_signal; zeros(Nsym/2,1)], sampsPerSym));
fltDelay = Nsym / (2*R); % Filter group delay, since raised cosine filter is linear phase and symmetric.�˲������ӳ٣���Ϊ�������˲�����������λ�ͶԳƵġ�
signal_complex = signal_complex(fltDelay*Fs+1:end); % Correct for propagation delay by removing filter transients ͨ��ȥ������˲�䣬��ȷ�ش����ӳ�
    
t = (0: length(signal_complex) - 1) / Fs;
f = linspace(-Fs/2,Fs/2,length(signal_complex));%�������Լ��ʸ��
I = real(signal_complex);
Q = imag(signal_complex);
FI = fftshift(fft(I));
FQ = fftshift(fft(Q));
a=I'.*cos(2*pi*fhz(i)*t);
b=Q'.*sin(2*pi*fhz(i)*t);
signal=signal+(a-b);
% signal= signal+(I'.*cos(2*pi*fhz(i)*t) - Q'.*sin(2*pi*fhz(i)*t));
end;

figure(4)
x=1:length(signal);
y1=abs(fft(signal,length(signal)));
plot(x(1:length(signal)/2),y1(1:length(signal)/2));
title('�����ź�Ƶ��');

figure(5)
x=1:length(signal);

plot(x,signal);
title('�����ź�ʱ����');
% figure(3);
% plot(f,abs(signal));

% Channel    
% noised_signal = awgn(signal,snr,'measured'); % Adding white gaussian noise
noised_signal = signal;
    
% Reciever
%����������OFDM�ź���ӵ��ź�


%ѡ���Ƿ�Ӵ�


%%���ݳ�ʼ��
ChannelNum = 50;        %MWCͨ����
N=40;
L=35;            %��ĿԤ��
M=35;
L0=floor(M/2);
TimeResolution=duration/length(noised_signal);
t_axis=0:TimeResolution:duration-TimeResolution;
fnyq=1/TimeResolution;
fp=fnyq/L;
Tp=1/fp;
m=ChannelNum;
K=81;
R=1;


%����+-1α�������
SignPatterns = randsrc(m,M);                


%%����ź�
MixedSigSequences = zeros(m,length(t_axis));
for channel=1:m
    MixedSigSequences(channel,:) = MixSignal(noised_signal,t_axis,SignPatterns(channel,:),Tp);
end

% figure(4)
% signal=MixedSigSequences(3,:)
% x=1:length(signal);
% y1=abs(fft(signal,length(signal)));
% plot(x(1:length(signal)/2),y1(1:length(signal)/2));
% title('����һ��֧·�ź�Ƶ��');

%% Analog low-pass filtering and actual sampling
% fprintf(1,'Filtering and decimation (=sampling)\n');
% ideal pass filter
temp = zeros(1,K);
temp(1) = 1;
lpf_z = interpft(temp,length(t_axis))/R/L; % impulse response

SignalSampleSequences = zeros(m,K);
NoiseSampleSequences = zeros(m,K);
% fprintf(1,'    Channel ');
decfactor = L*R;
for channel = 1:m
%     fprintf(1,'.');  if ( (mod(channel,5)==0) || (channel==m)) fprintf(1,'%d',channel); end
    SignalSequence =  MixedSigSequences(channel,:);
%     NoiseSequence   =  MixedNoiseSequences(channel,:);
    DigitalSignalSamples(channel, :) = FilterDecimate(SignalSequence,decfactor,lpf_z);
%     DigitalNoiseSamples(channel, :) = FilterDecimate(NoiseSequence,decfactor,lpf_z);
end
Digital_time_axis = downsample(t_axis,decfactor);
DigitalLength = length(Digital_time_axis);


%% CTF block
% fprintf(1,'---------------------------------------------------------------------------------------------\n');
% fprintf(1,'Entering CTF block\n');

% define matrices for fs=fp
S = SignPatterns;
theta = exp(-j*2*pi/L);
F = theta.^([0:L-1]'*[-L0:L0]);
np = 1:L0;
nn = (-L0):1:-1;
% This is for digital input only. Note that when R -> infinity,
% D then coincides with that of the paper
dn = [   (1-theta.^nn)./(1-theta.^(nn/R))/(L*R)      1/L    (1-theta.^np)./(1-theta.^(np/R))/(L*R)];
D = diag(dn);
A1 = S*F*D;
A = conj(A1);

DigitalSamples = DigitalSignalSamples;

% Frame construction
Q = DigitalSamples* DigitalSamples';
% decompose Q to find frame V
NumDomEigVals= FindNonZeroValues(eig(Q),5e-8);
[V,d] = eig_r(Q,min(NumDomEigVals,2*N));%V�Ƕ�Ӧ���������ֵ����λ�õ�����������d������2*N������ֵ
v = V*diag(sqrt(d));
% N iterations at most, since we force symmetry in the support...
[u, RecSupp] = RunOMP_Unnormalized(v, A, N, 0, 0.01, true);
RecSuppSorted = sort(unique(RecSupp));



%% Recover the singal
A_S = A(:,RecSuppSorted);
hat_zn = pinv(A_S)*DigitalSamples;  % inverting A_S
hat_zt = zeros(size(hat_zn,1),length(x));
for ii = 1:size(hat_zt,1)                     % interpolate (by sinc)
    hat_zt(ii,:) = interpft(hat_zn(ii,:),L*R*length(hat_zn(ii,:)));
end

x_rec = zeros(1,length(x));
for ii = 1:size(hat_zt,1)                      % modulate each band to their corresponding carriers
    x_rec = x_rec+hat_zt(ii,:).*exp(j*2*pi*(RecSuppSorted(ii)-L0-1)*fp.*t_axis);
end

x_rec = real(x_rec);

figure(5)
x=1:length(x_rec);
y1=abs(fft(x_rec,length(x_rec)));
plot(x(1:length(x_rec)/2),y1(1:length(x_rec)/2));
title('�ع������ź�Ƶ��');