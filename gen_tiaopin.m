function s=gen_tiaopin(fc,fs,N,len)

%% =============== ����FSK�ź� =====================
i=8;%�����ź���Ԫ��
n=16000;
a=round(rand(1,i));%�����������
t1=linspace(0,5e-3,n);
fc1=6e3;%�ز�1Ƶ��
fc2=8e3;%�ز�2Ƶ��
fm=i/5*1e3;%�����ź�Ƶ��
fc3=[];
b=[];
for ii = 1:8
    if(a(ii) >= 1)
        for p = 1:2000
            fc3 = [fc3 fc2];
            b = [b 1];
        end
    else
        for q = 1:2000
            fc3 = [fc3 fc1];
            b = [b 0];
        end
    end
end
fsk_sig = cos(2*pi*fc3.*t1);

L=fix(N/len);

s=[];
fc_sp=[1 25 59 93]*1e6;

for i=1:L
    k=unidrnd(length(fc_sp));
    ss=Generate_cw(fc+fc_sp(k),fs,len);
    s=[s ss];
end
% s = s.*fsk_sig;
if N-L*len>0
   for ii=1:N-L*len
    s=[s 0];
   end
end