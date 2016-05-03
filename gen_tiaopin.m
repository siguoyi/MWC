function [s ff]=gen_tiaopin(fc,fs,N,len)

%% =============== 生成FSK信号 =====================
i=8;%基带信号码元数
n=16000;
a=round(rand(1,i));%产生随机序列
t1=linspace(0,5e-3,n);
fc1=6e3;%载波1频率
fc2=8e3;%载波2频率
fc3=[];
ff=[];
b=[];
for ii = 1:8
    if(a(ii) >= 1)
        for p = 1:2000
            ff = [ff fc1];
            fc3 = [fc3 fc2];
            b = [b 1];
        end
    else
        for q = 1:2000
            ff = [ff fc1];
            fc3 = [fc3 fc1];
            b = [b 0];
        end
    end
end
fsk_sig = cos(2*pi*fc3.*t1);
L=fix(N/len);

s=[];
fc_sp=[11 25 59 93]*1e3;

one = ones(1,2000);
for i=1:L/2
    k=unidrnd(length(fc_sp));
    ss=Generate_cw(fc+fc_sp(k),fs,len);
    s=[s one ss];
end
s = s.*fsk_sig;
if N-L*len>0
   for ii=1:N-L*len
    s=[s 0];
    ff =[ff 0];
   end
end