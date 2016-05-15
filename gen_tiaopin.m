<<<<<<< HEAD
function [s fc1 fc2 s1 tt aa b]=gen_tiaopin(fc,fs,N,len)
=======
function [s fc1 fc2 s1 tt a]=gen_tiaopin(fc,fs,N,len)
>>>>>>> 197a28d402892d90384da5e3945ae482c1258284

%% =============== 生成FSK信号 =====================
num=64;;
i=num;%基带信号码元数
n=16000;
a=round(rand(1,i)) %产生随机序列
aa =a ;
t1=linspace(0,5e-3,n);
fc1=6e3;%载波1频率
fc2=8e3;%载波2频率
fc3=[];
tt=[t1];
s1=[];
b=[];
for ii = 1:num
    if(a(ii) >= 1)
        for p = 1:len
            fc3 = [fc3 fc2];
            b = [b 1];
        end
    else
        for q = 1:len
            fc3 = [fc3 fc1];
            b = [b 0];
        end
    end
end
fsk_sig = cos(2*pi*fc3.*t1);
% L=fix(N/len);
L=num;
s=[];
fc_sp=[11 25 59 93]*1e3;

% one = ones(1,2000);
for i=1:L
    k=unidrnd(length(fc_sp));
    ss=Generate_cw(fc+fc_sp(k),fs,len);
    s=[s ss];
end
s1 = [s];
s = s.*fsk_sig;
if N-L*len>0
   for ii=1:N-L*len
    s=[s 0];
    s1=[s1 0];
    tt=[tt 0];
   end
end