function  sig=generate_cw(fc,fs,N_signal)
t=[1:N_signal];
fnormal=fc/fs;
% sig=exp(j*2*pi*fnormal*t);
sig=cos(fnormal*t);