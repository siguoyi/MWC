%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo for Modulated Wideband Converter %
%             Version 1.0               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear,close all

%% Signal model
SNR = 10;                                   % Input SNR
N = 8;                                      % Number of bands (when counting  a band and its conjugate version separately)
B = 50e6;                                   % Maximal width of each band
Bi = ones(1,N/2)*B;
fnyq = 10e9;                                % Nyquist rate
Ei = rand(1,N/2)*10;                        % Energy of the i'th band
Tnyq = 1/fnyq;
R = 1;                                      % The length of the signal is R*(K+K0)*L
K = 91;
K0 = 10;                                    % R*K0*L is reserved for padding zeros
L = 195;
TimeResolution = Tnyq/R;
TimeWin = [0  L*R*K-1 L*R*(K+K0)-1]*TimeResolution; % Time interval in which signal is observed
Taui = [0.7 0.4 0.3]*max(TimeWin);          % Time offest of the i'th band

%% Sampling parameters
ChannelNum = 50;                            % Number of channels
L = 195;                                    % Aliasing rate
M = 195;
fp = fnyq/L;
fs = fp;                                    % Sampling rate at each channel, use fs=qfp, with odd q
m = ChannelNum;                             % Number of channels

% sign alternating  mixing
SignPatterns = randsrc(m,M);                % Draw a random +-1 for mixing sequences

% calucations
Tp = 1/fp;
Ts = 1/fs;
L0 = floor(M/2);                            
L = 2*L0+1;

%% Signal Representation
t_axis = TimeWin(1)  : TimeResolution : TimeWin(end);     % Time axis
t_axis_sig  = TimeWin(1)  : TimeResolution : TimeWin(2);

% Signal Generation
x = zeros(size(t_axis_sig));
fi = rand(1,N/2)*(fnyq/2-2*B) + B;      % Draw random carrier within [0, fnyq/2]
% for n=1:(N/2)
%     x = x+sqrt(Ei(n)) * sqrt(Bi(n))*sinc(Bi(n)*(t_axis_sig-Taui(n))) .* cos(2*pi*fi(n)*(t_axis_sig-Taui(n)));
% end
han_win = hann(length(x))';             % Add window
x = x.*han_win;
% x=real(exp(j*2*pi*10e6/100e6*([0:length(x)-1])));
[signal fc s1 tt] = gen_tiaopin(10e6,100e6,length(x),2000);
x=real(signal);
s1 = [s1 zeros(1,R*K0*L)];
plot(s1)
tt = [tt zeros(1,R*K0*L)];
x = [x, zeros(1,R*K0*L)];               % Zero padding
% figure(1)
% plot(x)
% figure(2)
% plot(abs(real(fft(x))))
%% Noise Generation
noise_nyq = randn(1,(K+K0)*L);              % Generate white Gaussian nosie within [-fnyq,fnyq]
noise = interpft(noise_nyq, R*(K+K0)*L);    % Interpolate into a finer grid (the same length as the signal)

% Calculate energies
NoiseEnergy = norm(noise)^2;
SignalEnergy = norm(x)^2;
CurrentSNR = SignalEnergy/NoiseEnergy;

%% Mixing
fprintf(1,'Mixing\n');

MixedSigSequences = zeros(m,length(t_axis));
for channel=1:m
    MixedSigSequences(channel,:) = MixSignal(x,t_axis,SignPatterns(channel,:),Tp);
end

MixedNoiseSequences = zeros(m,length(t_axis));
for channel=1:m
    MixedNoiseSequences(channel,:) = MixSignal(noise,t_axis,SignPatterns(channel,:),Tp);
end

%% Analog low-pass filtering and actual sampling
fprintf(1,'Filtering and decimation (=sampling)\n');
% ideal pass filter
temp = zeros(1,K+K0);
temp(1) = 1;
lpf_z = interpft(temp,length(t_axis))/R/L; % impulse response

SignalSampleSequences = zeros(m,K+K0);
NoiseSampleSequences = zeros(m,K+K0);
fprintf(1,'    Channel ');
decfactor = L*R;
for channel = 1:m
    fprintf(1,'.');  if ( (mod(channel,5)==0) || (channel==m)) fprintf(1,'%d',channel); end
    SignalSequence =  MixedSigSequences(channel,:);
    NoiseSequence   =  MixedNoiseSequences(channel,:);
    DigitalSignalSamples(channel, :) = FilterDecimate(SignalSequence,decfactor,lpf_z);
    DigitalNoiseSamples(channel, :) = FilterDecimate(NoiseSequence,decfactor,lpf_z);
end
Digital_time_axis = downsample(t_axis,decfactor);
DigitalLength = length(Digital_time_axis);

%% CTF block
fprintf(1,'---------------------------------------------------------------------------------------------\n');
fprintf(1,'Entering CTF block\n');

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
A = S*F*D;
A = conj(A);

SNR_val = 10^(SNR/10);          % not dB
% combine signal and noise
% DigitalSamples = DigitalSignalSamples + DigitalNoiseSamples*sqrt(CurrentSNR/SNR_val);
DigitalSamples = DigitalSignalSamples;

% Frame construction
Q = DigitalSamples* DigitalSamples';
% decompose Q to find frame V
NumDomEigVals= FindNonZeroValues(eig(Q),5e-8);
[V,d] = eig_r(Q,min(NumDomEigVals,2*N));
v = V*diag(sqrt(d));
% N iterations at most, since we force symmetry in the support...
% iterNum = size(Sorig,2);
[u, RecSupp] = RunOMP_Unnormalized(v, A, N, 0, 1e-12, false);
% [u, RecSupp] = RunSAMP_Unnormalized(v, A, 1);
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
% x_rec(:,1:2000) = 0;
% x_rec(:,4001:6000) = 0;
% x_rec(:,8001:10000) = 0;
% x_rec(:,12001:14000) = 0;
x_rec(:,16001:19695) = 0;
% sig = x + noise*sqrt(CurrentSNR/SNR_val);
snr1 = 20.*log10(norm(x(:,2601:3600))/norm(x(:,2601:3600)-x_rec(:,2601:3600)))
snr = 20.*log10(norm(x)/norm(x-x_rec))
%% demodulator
% fsk_sig = cos(2*pi*fc*tt(:,1:16000));
temp_fsk = x_rec(:,1:16000).*s1(:,1:16000);
[f,sf1] = T2F(tt(:,1:16000),temp_fsk);%通过低通滤波器
[t,st1] = lpf(f,sf1,1);
figure(1)
plot(st1)
fsk_sig = st1.*cos(2*pi*fc*tt(:,1:16000));
[f,sf2] = T2F(tt(:,1:16000),fsk_sig);%通过低通滤波器
[t,st2] = lpf(f,sf2,1);
figure(2)
plot(st2)
i=8;
for m=0:i-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%抽样判决
if st1(1,m*2000+1000)<5e-8;
for j=m*2000+1:(m+1)*2000;
at(1,j)=0;
end
else
for j=m*2000+1:(m+1)*20000;
at(1,j)=1;
end
end
end
figure(3)
plot(at)
%% plot module
% figure(1)
% subplot(211)
% plot(x)
% grid on;
% plot(x)
% title('Original signal');
% set(gca,'YLim',[-2 2]);
% subplot(212)
% plot(abs(real(fft(x))));
% title('Frequency spectrum of the original signal');
% 
% figure(2)
% subplot(211)
% plot(x)
% title('Original Signal');
% set(gca,'YLim',[-2 2]);
% subplot(212)
% plot(x_rec)
% title('Recostruction Signal');
% 
% figure(3)
% subplot(211)
% plot(abs(real(fft(x))));
% title('Frequency spectrum of original signal');
% subplot(212)
% plot(abs(real(fft(x_rec))));
% title('Frequency spectrum of recostruction signal');
% 
% figure(4)
% plot(abs(real(fft(x))));
% hold on
% plot(abs(real(fft(x_rec))),'r');
% legend('Original Signal', 'Recostruction Signal');
% title('Original Signal VS Recostruction Signal');
% 
% figure(5)
% subplot(311)
% plot(x(:,2001:4000))
% title('Original Signal');
% set(gca,'YLim',[-2 2]);
% subplot(312)
% plot(x_rec(:,2001:4000))
% title('Recostruction signal');
% subplot(313)
% plot(x(:,2001:4000))
% hold on;
% plot(x_rec(:,2001:4000),'r')