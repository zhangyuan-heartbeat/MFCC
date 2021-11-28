close all
clear
clc
[x]=audioread('FAML_Sa1.wav');
% MFCC implement with Matlab %  
bank=melbankm(24,256,16000,0,0.5,'t'); %Mel�˲����Ľ���Ϊ24��FFT�任�ĳ���Ϊ256������Ƶ��Ϊ16000Hz  
%��һ��Mel�˲�����ϵ��  
bank=full(bank); %full() convert sparse matrix to full matrix  
bank=bank/max(bank(:));  
for k=1:12  
    n=0:23;  
    dctcoef(k,:)=cos((2*n+1)*k*pi/(2*24));  
end  
w=1+6*sin(pi*[1:12]./12);%��һ��������������  
w=w/max(w);%Ԥ�����˲���  
xx=double(x);  
xx=filter([1-0.9375],1,xx);%�����źŷ�֡  
xx=enframe(xx,256,80);%��xx 256���Ϊһ֡  
%����ÿ֡��MFCC����  
for i=1:size(xx,1)  
    y=xx(i,:);  
    s=y'.*hamming(256);  
    t=abs(fft(s));%FFT���ٸ���Ҷ�任  
    t=t.^2;  
    c1=dctcoef*log(bank*t(1:129));  
    c2=c1.*w';  
    m(i,:)=c2;  
end  
%��һ�ײ��ϵ��  
dtm=zeros(size(m));  
for i=3:size(m,1)-2  
    dtm(i,:)=-2*m(i-2,:)-m(i-1,:)+m(i+1,:)+2*m(i+2,:);  
end  
dtm=dtm/3;  
%��ȡ���ײ��ϵ��  
dtmm=zeros(size(dtm));  
for i=3:size(dtm,1)-2  
    dtmm(i,:)=-2*dtm(i-2,:)-dtm(i-1,:)+dtm(i+1,:)+2*dtm(i+2,:);  
end  
dtmm=dtmm/3;  
%�ϲ�mfcc������һ�ײ��mfcc����  
ccc=[m dtm dtmm];  
%ȥ����β��֡����Ϊ����֡��һ�ײ�ֲ���Ϊ0  
ccc=ccc(3:size(m,1)-2,:);  
ccc; 
% subplot(2,1,1);  
% ccc_1=ccc(:,1); 
% m1=m(:,1)
figure(1)
plot(m);xlabel('Number of frames');ylabel('Amplitude'); 
figure(2)
plot(ccc);xlabel('Number of frames');ylabel('Amplitude');
  